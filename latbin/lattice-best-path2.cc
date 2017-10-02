// latbin/lattice-best-path2.cc

// Copyright (c) 2017 Joan Puigcerver <joapuipe@upv.es>

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <string>

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "fstext/fstext-utils2.h"
#include "fstext/lattice-weight.h"

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage =
        "Generate especial 1-best path through lattices, which minimizes a "
        "different loss-function other than the 0-1 error on the sequence.\n"
        "\n"
        "Traditional 1-best Viterbi decoding obtains:\n"
        "  argmax_y p(x, y) = argmax_y p(y | x)\n"
        "\n"
        "This minimizes the so called 0-1 loss at the sequence-level:\n"
        "  E_p(y|x) [ I(y != \\hat{y}) ]\n"
        "\n"
        "Instead, this decoding minimizes a different loss function which "
        "counts the errors on each position of the sequence, regardless of "
        "insertions and deletions:\n"
        "  E_p(y|x) [ \\sum_{k=1}^L I(y_k != \\hat{y}) + max(0, |y| - L) + "
        "max(0, |\\hat{y}| - L) ]\n"
        "\n"
        "where L = min(|y|, |\\hat{y}|).\n"
        "\n"
        "This loss is an upper bound on the expected Levenshtein distance, "
        "given p(y|x).\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Read(argc, argv);

    if (po.NumArgs() < 1 || po.NumArgs() > 2) {
      po.PrintUsage();
      exit(1);
    }

    const auto lattice_scale = LatticeScale(graph_scale, acoustic_scale);
    const std::string lattice_rspecifier = po.GetArg(1);
    Int32VectorWriter transcriptions_writer(po.GetOptArg(2));

    double total_cost = 0.0; int32 total_frames = 0;
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      // Read CompactLattice and disambiguate the input sequence length of
      // each state.
      // Note: I put this inside a different scope so the lattice can be freed
      // once it is modified.
      CompactLattice clat;
      std::vector<size_t> state_input_length;
      int32 lattice_frames = 0;
      {
        CompactLattice clat_tmp = lattice_reader.Value();
        lattice_reader.FreeCurrent();
        // Acoustic scale
        if (acoustic_scale != 1.0 || graph_scale != 1.0)
          ScaleLattice(lattice_scale, &clat_tmp);
        // Word insertion penalty
        if (insertion_penalty != 0.0)
          AddWordInsPenToCompactLattice(insertion_penalty, &clat_tmp);
        // If needed, sort the compact lattice in topological order
        TopSortCompactLatticeIfNeeded(&clat_tmp);
        // Compute state times. We don't need the vector, only the number of
        // frames in the lattice to compute the averages.
        std::vector<int32> state_times;
        lattice_frames = CompactLatticeStateTimes(clat_tmp, &state_times);
        KALDI_ASSERT(clat_tmp.Properties(kTopSorted, true) == kTopSorted);
        // Make sure that all sequences arriving to each state have the same
        // length. Note: Sort arcs in the lattice, to make sure that the
        // result is in topological order.
        ArcSort(&clat_tmp, OLabelCompare<CompactLatticeArc>());
        DisambiguateStateInputSequenceLength(
            clat_tmp, &clat, &state_input_length, false);
        KALDI_ASSERT(clat.Properties(kTopSorted, true) == kTopSorted);
        // Add disambiguation symbols at the end of the sequences so that
        // all transcriptions have the same length.
        AddSequenceLengthDismabiguationSymbol(&clat, &state_input_length);
        KALDI_ASSERT(clat.Properties(kTopSorted, true) == kTopSorted);
      }
      // Compute forward and backward LIKELIHOOD (- cost) of each state.
      std::vector<double> fw, bw;
      ComputeLatticeAlphasAndBetas(clat, false, &fw, &bw);

      typedef CompactLatticeArc::Label Label;
      typedef CompactLatticeArc::StateId StateId;
      std::map<std::tuple<Label, size_t>, double> word_pos_acc;
      for (StateIterator<CompactLattice> siter(clat); !siter.Done();
           siter.Next()) {
        const auto u = siter.Value();
        for (ArcIterator<CompactLattice> aiter(clat, u); !aiter.Done();
             aiter.Next()) {
          const auto& arc = aiter.Value();
          const auto v = arc.nextstate;
          const auto label = arc.olabel;
          if (label != 0) {
            const auto tup = std::make_tuple(label, state_input_length[v]);
            const double w = fw[u] + bw[v] - ConvertToCost(arc.weight);
            auto r = word_pos_acc.emplace(tup, w);
            if (!r.second) {
              r.first->second = LogAdd(r.first->second, w);
            }
          }
        }
      }

      // Normalize likelihood, so that we ket the posterior P(w | x, k) for
      // each position. The std::min() thing is to ignore precision errors that
      // make the posterior to be greater than 1.0.
      for (auto it = word_pos_acc.begin(); it != word_pos_acc.end(); ++it) {
        it->second = std::min(0.0, it->second - bw[clat.Start()]);
      }

      // Convert the lattice to a special fst in the tropical semiring,
      // where the weight in each arcs is equal to: 1 - P(w | x, k)
      VectorFst<StdArc> new_fst;
      for (StateId s = 0; s < clat.NumStates(); ++s) {
        new_fst.AddState();
      }
      new_fst.SetStart(clat.Start());

      for (StateIterator<CompactLattice> siter(clat); !siter.Done();
           siter.Next()) {
        const auto u = siter.Value();
        if (clat.Final(u) != CompactLatticeWeight::Zero()) {
          // Final states do not add any cost.
          new_fst.SetFinal(u, 0.0);
        }

        for (ArcIterator<CompactLattice> aiter(clat, u); !aiter.Done();
             aiter.Next()) {
          const auto& arc = aiter.Value();
          const auto v = arc.nextstate;
          const auto label = arc.olabel;
          if (label == 0) {
            // epsilon arcs do not add any cost
            new_fst.AddArc(u, StdArc(0, 0, 0.0, v));
          } else {
            const auto tup = std::make_tuple(label, state_input_length[v]);
            const double cost = exp(LogSub(0.0, word_pos_acc[tup]));
            new_fst.AddArc(u, StdArc(label, label, cost, v));
          }
        }
      }

      std::vector<int32> transcript_syms;
      double transcript_cost = 0.0;
      {
        VectorFst<StdArc> nbest_fst;
        ShortestPath(new_fst, &nbest_fst, 1);
        std::vector<VectorFst<StdArc>> nbests;
        ConvertNbestToVector(nbest_fst, &nbests);

        std::vector<StdArc::Label> isyms, osyms;
        StdArc::Weight cost;
        GetLinearSymbolSequence(nbests[0], &isyms, &osyms, &cost);
        transcript_cost = cost.Value();

        // Filter out epsilons and special disambiguation symbols.
        for (size_t k = 0; k < osyms.size(); ++k) {
          if (osyms[k] != 0 && osyms[k] != kNoLabel) {
            transcript_syms.push_back(osyms[k]);
          }
        }
      }

      if (po.NumArgs() > 1) {
        transcriptions_writer.Write(lattice_key, transcript_syms);
      }

      total_cost += transcript_cost;
      total_frames += lattice_frames;

      KALDI_LOG << "For utterance " << lattice_key << ", best cost is "
                << transcript_cost << " over " << lattice_frames << " frames.";
    }

    KALDI_LOG << "Overall cost per frame is " << (total_cost / total_frames)
              << " over " << total_frames << " frames.";

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

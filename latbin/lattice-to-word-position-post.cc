// latbin/lattice-to-word-position-post.cc

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

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage =
        "Compute the posterior log-probability of each word for each given "
        "transcription position. That is, we compute log P(w_k = v | x), "
        "for all possible transcript position k, and words v.\n"
        "You can optionally obtain the (best) frame segmentation of each "
        "word and position.\n"
        "\n"
        "Usage: lattice-to-word-position-post [options] lat-rspecifier "
        "post-wspecifier [segm-wspecifier]\n"
        " e.g.: lattice-to-word-position-post --acoustic-scale=0.1 ark:1.lats "
        "ark:1.word.pos.post\n"
        "See also: lattice-to-word-frame-post\n";

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

    if (po.NumArgs() != 2 || po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    const auto lattice_scale = LatticeScale(graph_scale, acoustic_scale);
    const std::string lattice_rspecifier = po.GetArg(1);
    PosteriorWriter posterior_writer(po.GetArg(2));
    // TODO: Write (best) segmentation for each (word, position)

    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();

      // Read CompactLattice and disambiguate the input sequence length of
      // each state.
      // Note: I put this inside a different scope so the lattice can be freed
      // once it is modified.
      CompactLattice clat;
      std::vector<size_t> state_input_length;
      size_t max_len = 0;
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

        // Make sure that all sequences arriving to each state have the same
        // length.
        max_len = DisambiguateStateInputSequenceLength(
            clat_tmp, &clat, &state_input_length, false);
      }
      // Compute forward and backward likelihoods of each state in the
      // (modified) lattice.
      std::vector<double> fw_lkh, bw_lkh;
      const double total_lkh =
          ComputeLatticeAlphasAndBetas(clat, false, &fw_lkh, &bw_lkh);
      // Compute total likelihood for each (word, position).
      // Note: acc[0] is not used.
      std::vector<std::map<CompactLatticeArc::Label, double>> acc(max_len + 1);
      for (StateIterator<CompactLattice> siter(clat); !siter.Done();
           siter.Next()) {
        const auto u = siter.Value();
        for (ArcIterator<CompactLattice> aiter(clat, u);
             !aiter.Done(); aiter.Next()) {
          const auto& arc  = aiter.Value();
          const auto& word = arc.ilabel;
          const size_t pos = state_input_length[arc.nextstate];
          if (word != 0) {
            const double total_lkh_through_arc =
                fw_lkh[u] + bw_lkh[arc.nextstate]
                -(arc.weight.Weight().Value1() + arc.weight.Weight().Value2());
            auto r = acc[pos].emplace(word, total_lkh_through_arc);
            if (!r.second) {
              r.first->second = LogAdd(r.first->second, total_lkh_through_arc);
            }
          }
        }
      }
      // Normalize likelihoods to get posteriors.
      PosteriorWriter::T posterior(max_len);
      for (size_t n = 1; n <= max_len; ++n) {
        for (const auto& kv : acc[n]) {
          const auto label = kv.first;
          const auto logp  = kv.second - total_lkh;
          posterior[n - 1].emplace_back(label, logp);
        }
        // Sort posteriors within a position in decreasing order.
        std::sort(posterior[n - 1].begin(), posterior[n - 1].end(),
                  [](const std::pair<int32, BaseFloat>& a,
                     const std::pair<int32, BaseFloat>& b) -> bool {
                    if (a.second != b.second) return b.second < a.second;
                    else return a.first < b.first;
                  });
      }
      // Write posteriors.
      posterior_writer.Write(lattice_key, posterior);
    }
    posterior_writer.Close();
    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

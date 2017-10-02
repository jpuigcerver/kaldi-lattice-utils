// kwsbin2/lattice-index.cc

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
#include "util/basic-tuple-vector-holder.h"

#include "fst/script/print.h"

namespace fst {

template <typename Arc>
void CreateQueryFst(const std::vector<typename Arc::Label>& alphabet,
                    const typename Arc::Label& query_label,
                    MutableFst<Arc>* fst) {
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;

  fst->DeleteStates();
  fst->AddState();
  fst->AddState();
  fst->SetStart(0);

  for (const Label& label : alphabet) {
    if (label == 0) continue;
    if (label == query_label) {
      fst->AddArc(0, Arc(label, label, Weight::One(), 1));
    } else {
      fst->AddArc(0, Arc(label, label, Weight::One(), 0));
    }
    fst->AddArc(1, Arc(label, label, Weight::One(), 1));
  }

  fst->SetFinal(0, Weight::Zero());
  fst->SetFinal(1, Weight::One());

  constexpr uint64_t properties =
      kExpanded | kMutable | kAcceptor | kIDeterministic | kODeterministic |
      kNoEpsilons | kNoIEpsilons | kNoOEpsilons | kUnweighted |
      kInitialCyclic | kAccessible | kCoAccessible;
  fst->SetProperties(kFstProperties, properties);
}

}  // namespace fst

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage =
        "Create an inverted index of the given lattices, where the score of "
        "each word is the probability that the word is in the transcription "
        "of the signal. I.e. P(R = 1 | x, v), where x is the image and v is "
        "the word.\n"
        "\n"
        "For additional details, check the paper \"Probabilistic "
        "interpretation and improvements to the HMM-filler for handwritten "
        "keyword spotting\", by J. Puigcerver et al.\n"
        "\n"
        "Usage: lattice-index [options] <lattice-rspecifier> "
        "<index-wspecifier>\n"
        " e.g.: lattice-index --acoustic-scale=0.1 ark:1.lats ark:1.index\n";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
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
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Read(argc, argv);

    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    typedef TableWriter<BasicTupleVectorHolder<int32, double>>
        UtteranceScoreWriter;

    const auto lattice_scale = LatticeScale(graph_scale, acoustic_scale);
    const std::string lattice_rspecifier = po.GetArg(1);
    UtteranceScoreWriter utterance_writer(po.GetArg(2));

    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(lattice_scale, &clat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(insertion_penalty, &clat);
      // Lattice prunning
      if (beam != std::numeric_limits<BaseFloat>::infinity())
        PruneLattice(beam, &clat);
      // If needed, sort the compact lattice in topological order
      TopSortCompactLatticeIfNeeded(&clat);

      // Compute lattice total likelihood
      std::vector<double> bw_lkh;
      ComputeCompactLatticeBetas(clat, &bw_lkh);
      const double total_lkh = bw_lkh[clat.Start()];
      // Obtain the list of words in the hypotheses set.
      std::vector<int32> symbols;
      GetOutputSymbols(clat, false /* do not include epsilon */, &symbols);

      std::vector<std::tuple<int32, double>> utterance_scores;
      for (const auto& symbol : symbols) {
        CompactLattice query_lat, composed_lat;
        // Create the fst representing the language \Sigma^* symbol \Sigma^*
        // The query fst must be deterministic so that no repeated paths are
        // added in the composition with the log_fst (i.e. original lattice).
        CreateQueryFst(symbols, symbol, &query_lat);
        // Compose / intersect with the lattice with the query fst to obtain
        // all paths in the lattice containing the query word.
        Compose(clat, query_lat, &composed_lat);
        TopSortCompactLatticeIfNeeded(&composed_lat);
        ComputeCompactLatticeBetas(composed_lat, &bw_lkh);
        const double query_lkh = bw_lkh[composed_lat.Start()];
        if (query_lkh > total_lkh) {
          KALDI_WARN << "Found query with higher likelihood than the "
                     << "likelihood (" << query_lkh << " vs. "
                     << total_lkh << ")";
        }
        utterance_scores.emplace_back(symbol, query_lkh - total_lkh);
      }
      // Sort scores in decreasing order.
      std::sort(utterance_scores.begin(), utterance_scores.end(),
                [](const std::tuple<int32, double>& a,
                   const std::tuple<int32, double>& b) -> bool {
                  if (std::get<1>(b) != std::get<1>(a)) {
                    return std::get<1>(b) < std::get<1>(a);
                  } else {
                    return std::get<0>(a) < std::get<0>(b);
                  }
                });
      // Write scores to the table
      utterance_writer.Write(lattice_key, utterance_scores);
    }
    utterance_writer.Close();
    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

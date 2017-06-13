// latbin/lattice-to-transcript-length-post.cc

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
        "Compute the distribution of the length of the transcriptions in a "
        "lattice.\n"
        "\n"
        "Usage: lattice-to-transcript-length-post [options] "
        "lattice-rspecifier1 posterior-wspecifier\n"
        " e.g.: lattice-to-transcript-length-post ark:1.lats ark:1.posts\n";

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

    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    const auto lattice_scale = LatticeScale(graph_scale, acoustic_scale);
    const std::string lattice_rspecifier = po.GetArg(1);
    PosteriorWriter posterior_writer(po.GetArg(2));

    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();

      // Read CompactLattice and disambiguate the input sequence length of
      // each state.
      // Note: I put this inside a different scope so the lattice can be freed
      // once it is modified.
      CompactLattice clat;
      std::vector<size_t> state_input_length;
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
        DisambiguateStateInputSequenceLength(
            clat_tmp, &clat, false, &state_input_length);
      }
      // Compute forward and backward likelihoods of each state in the
      // (modified) lattice.
      std::vector<double> fw_lkh, bw_lkh;
      const double total_lkh =
          ComputeLatticeAlphasAndBetas(clat, false, &fw_lkh, &bw_lkh);
      // Accumulate the likelihoods of ending in each state of the same length,
      // in order to compute the total likelihood of each transcript length.
      std::map<size_t, double> acc;
      for (StateIterator<CompactLattice> siter(clat); !siter.Done();
           siter.Next()) {
        const auto u = siter.Value();
        if (clat.Final(u) != CompactLatticeWeight::Zero()) {
          const size_t length = state_input_length[u];
          const double lkh_state_end = fw_lkh[u] - ConvertToCost(clat.Final(u));
          auto r = acc.emplace(length, lkh_state_end);
          if (!r.second) {
            r.first->second = LogAdd(r.first->second, lkh_state_end);
          }
        }
      }
      // Normalize likelihoods to get posteriors.
      PosteriorWriter::T posterior(1);
      for (const auto& kv : acc) {
        const auto length = kv.first;
        const auto logp   = kv.second - total_lkh;
        posterior[0].emplace_back(length, logp);
      }
      // Sort posteriors within a position in decreasing order.
      std::sort(posterior[0].begin(), posterior[0].end(),
                [](const std::pair<int32, BaseFloat>& a,
                   const std::pair<int32, BaseFloat>& b) -> bool {
                  if (a.second != b.second) return b.second < a.second;
                  else return a.first < b.first;
                });
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

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
#include "base/timer.h"
#include "fstext/kaldi-fst-io.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/common-utils.h"


int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage =
        "Compute the posterior log-probability of each word for each given "
        "utterance frame. That is, we compute log P(a_i = v | x), for all "
        "possible utterance frames i and words v.\n"
        "\n"
        "Usage: lattice-to-word-frame-post [options] lat-rspecifier "
        "post-wspecifier\n"
        " e.g.: lattice-to-word-frame-post --acoustic-scale=0.1 ark:1.lats "
        "ark:1.word.pos.post\n"
        "See also: lattice-to-word-position-post\n";

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
    size_t num_lattices = 0;
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next(), ++num_lattices) {
      const std::string lattice_key = lattice_reader.Key();
      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(lattice_scale, &clat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(insertion_penalty, &clat);
      // If needed, sort the compact lattice in topological order
      TopSortCompactLatticeIfNeeded(&clat);
      Timer timer;
      // Compute the times of each state.
      // Assumption: The lattice must be aligned!
      std::vector<int32> times;
      const int32 total_frames = kaldi::CompactLatticeStateTimes(clat, &times);
      // Compute forward and backward likelihoods of each state in the
      // (modified) lattice.
      std::vector<double> fw_lkh, bw_lkh;
      const double total_lkh =
          ComputeLatticeAlphasAndBetas(clat, false, &fw_lkh, &bw_lkh);
      // Compute total likelihood for each (word, frame).
      bool misaligned_lattice = false;
      std::vector<std::map<CompactLatticeArc::Label, double>> acc(total_frames);
      for (StateIterator<CompactLattice> siter(clat); !siter.Done();
           siter.Next()) {
        const auto u = siter.Value();
        for (ArcIterator<CompactLattice> aiter(clat, u);
             !aiter.Done(); aiter.Next()) {
          const auto& arc = aiter.Value();
          const auto word = arc.ilabel;
          if (word != 0) {
            const double total_lkh_through_arc =
                fw_lkh[u] + bw_lkh[arc.nextstate]
                -(arc.weight.Weight().Value1() + arc.weight.Weight().Value2());
            if (times[u] >= times[arc.nextstate] && !misaligned_lattice) {
              misaligned_lattice = true;
              KALDI_WARN << "Lattice " << lattice_key << " is misaligned, "
                         << "a word with zero duration was found!";
            }
            for (int32 k = times[u]; k < times[arc.nextstate]; ++k) {
              auto r = acc[k].emplace(word, total_lkh_through_arc);
              if (!r.second) {
                r.first->second =
                    LogAdd(r.first->second, total_lkh_through_arc);
              }
            }
          }
        }
      }
      // Normalize likelihoods to get posteriors.
      PosteriorWriter::T posterior(total_frames);
      for (size_t n = 0; n < total_frames; ++n) {
        for (const auto& kv : acc[n]) {
          const auto label = kv.first;
          const auto logp  = kv.second - total_lkh;
          posterior[n].emplace_back(label, logp);
        }
        // Sort posteriors within a frame in decreasing order.
        std::sort(posterior[n].begin(), posterior[n].end(),
                  [](const std::pair<int32, BaseFloat>& a,
                     const std::pair<int32, BaseFloat>& b) -> bool {
                    if (a.second != b.second) return b.second < a.second;
                    else return a.first < b.first;
                  });
      }
      KALDI_LOG << "Lattice " << lattice_key << " posteriorgram computed in "
                << timer.Elapsed() << " seconds.";
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

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

int main(int argc, char** argv) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage = "";

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
      // Make sure that all sequences arriving to each state have the same
      // length.
      max_len = DisambiguateStateInputSequenceLength(
          clat_tmp, &clat, &state_input_length, false);
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
}

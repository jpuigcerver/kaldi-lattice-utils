// latbin/lattice-restrict-length.cc

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
        "Restrict paths from lattices such that the length of the "
        "transcription is the one specified (read from a separate table).\n"
        "\n"
        "Usage: lattice-restrict-length [options] <length-rspecifier> "
        "<lattice-rspecifier> <lattice-wspecifier>\n"
        " e.g.: lattice-restrict-length ark:1.lengths ark:1.lats ark:1r.lats\n";

    ParseOptions po(usage);
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    RandomAccessTableReader<BasicHolder<uint32_t>> length_reader(po.GetArg(1));
    SequentialCompactLatticeReader lattice_reader(po.GetArg(2));
    CompactLatticeWriter lattice_writer(po.GetArg(3));

    size_t n_done = 0, n_fail = 0;
    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      if (!length_reader.HasKey(lattice_key)) {
        KALDI_WARN << "Lattice key " << lattice_key
                   << " was not found in the length table archive.";
        n_fail++;
        continue;
      }

      // Read CompactLattice and disambiguate the input sequence length of
      // each state. Meaning that all sequences arriving to each state have
      // the same length.
      CompactLattice clat;
      std::vector<size_t> state_input_length;
      {
        CompactLattice clat_tmp = lattice_reader.Value();
        lattice_reader.FreeCurrent();
        // If needed, sort the compact lattice in topological order
        TopSortCompactLatticeIfNeeded(&clat_tmp);
        // Make all sequences arriving to each state have the same length.
        DisambiguateStateInputSequenceLength(
            clat_tmp, &clat, &state_input_length, false);
      }
      const size_t given_length = length_reader.Value(lattice_key);
      // Since we are given the true length of the transcription, set the final
      // weight of all states with input_length != given_length to zero.
      for (StateIterator<CompactLattice> siter(clat); !siter.Done();
           siter.Next()) {
        const auto u = siter.Value();
        if (state_input_length[u] != given_length) {
          clat.SetFinal(u, CompactLatticeWeight::Zero());
        }
      }
      // Trim the lattice removing useless arcs and states.
      Connect(&clat);

      // Write
      lattice_writer.Write(lattice_key, clat);
      n_done++;
    }

    KALDI_LOG << "Done " << n_done << " lattices, failed for " << n_fail;
    return (n_done != 0);
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

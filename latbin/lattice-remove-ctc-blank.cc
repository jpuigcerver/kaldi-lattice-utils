// MIT License
//
// Copyright (c) 2016 Joan Puigcerver <joapuipe@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "fstext/determinize-lattice.h"

namespace kaldi {

void RemoveCTCBlankFromLattice(
    const Lattice& inp, const LatticeArc::Label blank, Lattice* out) {
  typedef LatticeArc::StateId StateId;
  typedef LatticeArc::Label Label;
  // Get table mapping from output symbols from the output labels in the input
  // lattice to states
  std::unordered_map<Label, StateId> symbol2state;
  for (fst::StateIterator<Lattice> siter(inp); !siter.Done(); siter.Next()) {
    const StateId s = siter.Value();
    for (fst::ArcIterator<Lattice> aiter(inp, s); !aiter.Done();
         aiter.Next()) {
      const LatticeArc& arc = aiter.Value();
      const Label o = arc.olabel;
      if (o != blank && o != 0 && symbol2state.count(o) == 0) {
        symbol2state.insert(std::make_pair(o, symbol2state.size() + 1));
      }
    }
  }
  // Create composition lattice, such that output = compose(input, C)
  Lattice C;
  for (size_t s = 0; s < symbol2state.size() + 1; ++s) {
    C.SetFinal(C.AddState(), LatticeWeight::One());
  }
  C.SetStart(0);
  // Self-loop in the blank state
  C.AddArc(0, LatticeArc(blank, 0, LatticeWeight::One(), 0));
  for (const std::pair<Label, StateId>& p : symbol2state) {
    // Arc from the initial state to the symbol's state
    C.AddArc(0, LatticeArc(p.first, p.first, LatticeWeight::One(), p.second));
    // Self-loop in the symbol's state
    C.AddArc(p.second, LatticeArc(p.first, 0, LatticeWeight::One(), p.second));
    // Arc back to the blank state
    C.AddArc(p.second, LatticeArc(blank, 0, LatticeWeight::One(), 0));
    // Arc to all other symbol states
    for (const std::pair<Label, StateId>& p2 : symbol2state) {
      if (p.first != p2.first) {
        C.AddArc(p.second, LatticeArc(p2.first, p2.first,
                                      LatticeWeight::One(), p2.second));
      }
    }
  }
  // Compute output lattice
  fst::Compose(inp, C, out);
}

void ProcessLattice(const std::string& key, Lattice *ilat, Lattice *olat,
                    const LatticeArc::Label blank_symbol, const double beam,
                    const double acoustic_scale, const double graph_scale,
                    const bool only_best_alignment) {
  // Make sure that lattice complies with all assumptions
  const uint64_t properties = ilat->Properties(fst::kAcceptor | fst::kAcyclic,
                                               true);
  if ((properties & fst::kAcceptor) != fst::kAcceptor) {
    KALDI_ERR << "Lattice " << key << " is not an acceptor";
  }
  if ((properties & fst::kAcyclic) != fst::kAcyclic) {
    KALDI_ERR << "Lattice " << key << " is not acyclic";
  }

  // Prune lattice
  if (beam != std::numeric_limits<BaseFloat>::infinity()) {
    // Scale lattice
    if (acoustic_scale != 1.0 || graph_scale != 1.0) {
      std::vector<std::vector<double> > scale(2, std::vector<double>(2, 0));
      scale[0][0] = graph_scale;
      scale[1][1] = acoustic_scale;
      fst::ScaleLattice(scale, ilat);
    }

    PruneLattice(beam, ilat);

    // Put lattices in the original scale
    if (acoustic_scale != 1.0 || graph_scale != 1.0) {
      std::vector<std::vector<double> > inv_scale(2, std::vector<double>(2, 0));
      inv_scale[0][0] = 1.0 / graph_scale;
      inv_scale[1][1] = 1.0 / acoustic_scale;
      fst::ScaleLattice(inv_scale, ilat);
    }
  }

  // Remove CTC Blanks from the output symbols
  Lattice tmp;
  RemoveCTCBlankFromLattice(*ilat, blank_symbol,
                            only_best_alignment ? &tmp : olat);

  // Determinize to keep only the best alignment hypothesis for each sequence
  // of characters.
  if (only_best_alignment) {
    fst::Invert(&tmp);
    fst::DeterminizeLattice<LatticeWeight, int32>(tmp, olat);
    fst::Invert(olat);
  }
}

}  // namespace kaldi

int main(int argc, char** argv) {
  try {
    using namespace kaldi;

    const char* usage =
        "Remove CTC blank symbols from the output labels of Kaldi lattices.\n"
        "\n"
        "Usage: lattice-remove-ctc-blank blank-symbol lat-rspecifier lat-wspecifier\n"
        " e.g.: lattice-remove-ctc-blank 32 ark:input.ark ark:output.ark\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    bool only_best_alignment = false;
    bool write_compact = true;
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                               "and adding the insertion penalty).");
    po.Register("only-best-alignment", &only_best_alignment,
                "If true, keep only the most likely alignment for each "
                "sequence of characters.");
    po.Register("write-compact", &write_compact,
                "If true, write compact lattices.");
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    const std::string blank_symbol_str = po.GetArg(1);
    const std::string lattice_in_str = po.GetArg(2);
    const std::string lattice_out_str = po.GetArg(3);
    const bool lattice_in_is_table =
        (ClassifyRspecifier(lattice_in_str, NULL, NULL) != kNoRspecifier);
    const bool lattice_out_is_table =
        (ClassifyWspecifier(lattice_out_str, NULL, NULL, NULL) != kNoWspecifier);

    LatticeArc::Label blank_symbol = 0;
    if (!ConvertStringToInteger(blank_symbol_str, &blank_symbol)) {
      KALDI_ERR << "String \"" << blank_symbol_str
                << "\" cannot be converted to an integer";
    }
    if (blank_symbol == 0) {
      KALDI_ERR << "Symbol 0 is reserved for epsilon!";
    }


    if (lattice_in_is_table && lattice_out_is_table) {
      SequentialLatticeReader lattice_reader(lattice_in_str);

      CompactLatticeWriter *clattice_writer = nullptr;
      LatticeWriter *lattice_writer = nullptr;
      if (write_compact) {
        clattice_writer = new CompactLatticeWriter(lattice_out_str);
      } else {
        lattice_writer = new LatticeWriter(lattice_out_str);
      }
      for (; !lattice_reader.Done(); lattice_reader.Next()) {
        // Read input lattice
        const std::string lattice_key = lattice_reader.Key();
        Lattice ilat = lattice_reader.Value(), olat;
        lattice_reader.FreeCurrent();
        ProcessLattice(lattice_key, &ilat, &olat, blank_symbol, beam,
                       acoustic_scale, graph_scale, only_best_alignment);
        if (write_compact) {
          Lattice tmp;
          fst::Push<Lattice::Arc, fst::REWEIGHT_TO_INITIAL>(olat, &tmp, fst::kPushLabels);
          CompactLattice colat;
          ConvertLattice(tmp, &colat, /*invert=*/true);
          clattice_writer->Write(lattice_key, colat);
          KALDI_LOG << olat.NumStates() << " " << colat.NumStates();
        } else {
          lattice_writer->Write(lattice_key, olat);
        }
      }
      if (write_compact) {
        delete clattice_writer;
      } else {
        delete lattice_writer;
      }
    } else {
      KALDI_ERR << "Not implemented! Both input and output lattices must be "
                << "Kaldi tables.";
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}

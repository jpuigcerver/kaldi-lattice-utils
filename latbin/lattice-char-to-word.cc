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

#include "fstext/expand_subpaths_fst.h"

namespace fst {

template<typename M>
void MapToSymbolsTable(
    const M &smap, SymbolTable *stable) {
  typedef typename M::key_type key_type;        // i.e. vector<int>
  typedef typename M::mapped_type mapped_type;  // i.e. int
  // First, put map into a vector, in order to sort it according to the values,
  // not the keys.
  std::vector <std::pair<mapped_type, std::string>> vmap;
  for (const auto &k_v : smap) {
    const key_type &k = k_v.first;
    const mapped_type &v = k_v.second;
    std::ostringstream ss;
    for (const mapped_type &c : k) {
      ss << c << "_";
    }
    std::string str = ss.str();
    if (str.empty()) { str = "0"; }  // epsilon
    else { str.pop_back(); }         // pop last _
    vmap.push_back(std::make_pair(v, str));
  }
  // Sort the vector.
  std::sort(vmap.begin(), vmap.end());
  // Add the symbols to the output table.
  for (const auto &k_v : vmap) {
    const mapped_type &k = k_v.first;
    const std::string &v = k_v.second;
    const mapped_type &kt = stable->AddSymbol(v, k);
    KALDI_ASSERT(kt == k);
  }
}

}  // namespace fst



int main(int argc, char **argv) {
  try {
    using namespace kaldi;
    typedef CompactLattice::Label Label;

    const char *usage =
        "Convert character-level lattices into word-level lattices by "
            "expanding the subpaths in between any of two separator symbols.\n"
            "\n"
            "Keep in mind that this lattice expansion has an exponential cost. "
            "For instance, if the set of separator symbols was empty, all paths "
            "from the input lattice would be expanded, so that each arc in the "
            "output lattice would be a full path from the input lattice.\n"
            "\n"
            "However, the exponential growth is constrained by the use of "
            "separator symbols, and make the tool practical in real scenarios.\n"
            "\n"
            "In addition, there are two pruning mechanisms to prevent the output "
            "lattices from exploding:\n"
            "\n"
            "1. You can prune the character lattices before expanding with the "
            "--beam option.\n"
            "2. You can set a maximum length for the output words with the "
            "--max-length option. Any path with a word longer than this number "
            "of characters will be removed from the output path.\n"
            "\n"
            "Usage: lattice-char-to-word [options] separator-symbols "
            "lat-rspecifier lat-wspecifier\n"
            " e.g.: lattice-char-to-word \"3 4\" ark:1.lat ark:1-words.lat\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    int32 max_length = std::numeric_limits<int32>::max();
    std::string save_symbols = "";
    std::string other_groups_str = "";

    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Register("save-symbols", &save_symbols, "If given, all lattices will "
        "use the same symbol table which will be written to this "
        "destination file. If not provided, each lattice contains "
        "its own symbol table.");
    po.Register("max-length", &max_length,
                "Max. length (in characters) for a word.");
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    if (graph_scale <= 0.0 || acoustic_scale <= 0.0) {
      KALDI_ERR
          << "--acoustic-scale and --graph-scale must be strictly greater than 0.0!";
    }


    LabelGroup<Label> label_group;
    label_group.AddGroupsStr(po.GetArg(1));
    const std::string lattice_in_str = po.GetArg(2);
    const std::string lattice_out_str = po.GetArg(3);

    std::unordered_map<std::vector<Label>, Label> label_map;
    SequentialCompactLatticeReader lattice_reader(lattice_in_str);
    CompactLatticeWriter lattice_writer(lattice_out_str);

    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();

      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();

      // Prune lattice
      if (beam != std::numeric_limits<BaseFloat>::infinity()) {
        // Acoustic scale
        if (acoustic_scale != 1.0 || graph_scale != 1.0)
          ScaleLattice(fst::LatticeScale(graph_scale, acoustic_scale), &clat);
        PruneLattice(beam, &clat);
        // Put weights into the original scale
        if (acoustic_scale != 1.0 || graph_scale != 1.0)
          ScaleLattice(fst::LatticeScale(1.0 / graph_scale,
                                         1.0 / acoustic_scale), &clat);
      }

      // Expand fst to obtain the word-level arcs
      CompactLattice olat;
      if (save_symbols == "") label_map.clear();

      ExpandSubpathsByLabelGroup(clat, label_group, &olat, &ilabel_map, &olabel_map);


      // Write symbol table in the output lattice, if no filename was given
      // for a common table.
      if (save_symbols == "") {
        fst::SymbolTable stable;
        MapToSymbolsTable(label_map, &stable);
        olat.SetInputSymbols(&stable);
        olat.SetOutputSymbols(&stable);
      }

      // Write output lattice
      lattice_writer.Write(lattice_key, olat);
    }

    // Save symbol table into separate files, if these
    // filenames where given
    if (save_symbols != "") {
      fst::SymbolTable stable;
      MapToSymbolsTable(label_map, &stable);
      const bool r = stable.WriteText(save_symbols);
      KALDI_ASSERT(r);
    }

    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return 1;
  }
}

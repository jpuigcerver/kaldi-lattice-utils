// MIT License
//
// Copyright (c) 2016-2018 Joan Puigcerver <joapuipe@gmail.com>
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
#include "base/timer.h"
#include "fstext/expand-subpaths-same-label-class.h"
#include "fstext/fstext-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/label-group.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"

namespace kaldi {

using fst::ExpandSubpathsOptions;
using fst::ExpandSubpathsWithSameLabelClass;
using fst::LabelGroup;
using fst::SymbolTable;

template<typename L, typename LW>
class ExpandLatticeTask {
 public:
  typedef typename L::Arc Arc;
  typedef typename L::Arc::Label Label;

  ExpandLatticeTask(const std::string &key,
                    const L &lat,
                    const LabelGroup<Label> &label_group,
                    const BaseFloat &acoustic_scale,
                    const BaseFloat &graph_scale,
                    const BaseFloat &beam,
                    const ExpandSubpathsOptions &expand_opts,
                    LW *lattice_writer,
                    SymbolTable *symbol_table)
      : key_(key),
        lat_(lat),
        label_group_(label_group),
        acoustic_scale_(acoustic_scale),
        graph_scale_(graph_scale),
        beam_(beam),
        opts_(expand_opts),
        lattice_writer_(lattice_writer),
        symbol_table_(symbol_table) {}

  ~ExpandLatticeTask() {
    // Note: Destructor is executed sequentially
    // Copy symbols to the given symbol table
    if (symbol_table_) {
      if (lat_.InputSymbols()) {
        symbol_table_->AddTable(*lat_.InputSymbols());
      }
      if (lat_.OutputSymbols()) {
        symbol_table_->AddTable(*lat_.OutputSymbols());
      }
      fst::Relabel(&lat_, symbol_table_, symbol_table_);
      lat_.SetInputSymbols(nullptr);
      lat_.SetOutputSymbols(nullptr);
    }
    // Write output lattice
    lattice_writer_->Write(key_, lat_);
  }

  void operator()() {
    Timer timer;

    auto orig_num_states = lat_.NumStates();
    auto orig_num_arcs = NumArcs(lat_);

    // Prune lattice
    if (beam_ != std::numeric_limits<BaseFloat>::infinity()) {
      // Acoustic scale
      if (acoustic_scale_ != 1.0 || graph_scale_ != 1.0)
        fst::ScaleLattice(fst::LatticeScale(graph_scale_,
                                            acoustic_scale_), &lat_);
      PruneLattice(beam_, &lat_);
      // Put weights into the original scale
      if (acoustic_scale_ != 1.0 || graph_scale_ != 1.0)
        fst::ScaleLattice(fst::LatticeScale(1.0 / graph_scale_,
                                            1.0 / acoustic_scale_), &lat_);

      const auto aux = NumArcs(lat_);
      KALDI_VLOG(1) << "Lattice " << key_
                    << " pruned #states from " << orig_num_states
                    << " to " << lat_.NumStates()
                    << " and #arcs from " << orig_num_arcs
                    << " to " << aux;
      orig_num_arcs = aux;
      orig_num_states = lat_.NumStates();
    }

    // Note: Optionally, report the time spend in the preparing the lattice
    // alone. This time is included in the total cost of the lattice.
    KALDI_VLOG(1) << "Lattice " << key_ << " preprocessed in "
                  << timer.Elapsed() << " seconds.";

    // Expansion of subpaths formed by arcs of the same group of labels.
    // WARNING: THIS HAS AN EXPONENTIAL WORST CASE COST
    const auto out_num_arcs =
        ExpandSubpathsWithSameLabelClass<Arc, Label>(label_group_, &lat_, opts_);
    const auto out_num_states = lat_.NumStates();
    KALDI_LOG << "Lattice " << key_ << " expanded #states from "
              << orig_num_states << " to " << out_num_states
              << " and #arcs from " << orig_num_arcs << " to " << out_num_arcs
              << " in " << timer.Elapsed() << " seconds.";
    KALDI_VLOG(1) << "Lattice " << key_
                  << ": input subpath symbols = "
                  << lat_.InputSymbols()->NumSymbols()
                  << ", output subpath symbols = "
                  << lat_.OutputSymbols()->NumSymbols();
  }

 protected:
  const std::string key_;
  L lat_;
  LabelGroup<Label> label_group_;
  BaseFloat acoustic_scale_, graph_scale_, beam_;
  ExpandSubpathsOptions opts_;
  LW *lattice_writer_;
  SymbolTable *symbol_table_;
};

}  // namespace kaldi

namespace fst {

bool WriteSymbolTable(const SymbolTable &table,
                      const std::string& filename,
                      bool text) {
  if (text) {
    return table.WriteText(filename);
  } else {
    return table.Write(filename);
  }
}

SymbolTable* ReadOrCreateSymbolTable(const std::string& filename) {
  SymbolTable* table = SymbolTable::Read(filename);
  if (!table) {
    table = new SymbolTable;
  }
  return table;
}

}  // namespace fst


int main(int argc, char **argv) {
  try {
    using namespace kaldi;
    typedef CompactLattice::Arc::Label Label;

    const char *usage =
        "Expand subpaths in lattices where all labels in the path have "
        "the same \"class\".\n"
        "\n"
        "This is useful, for instance, to convert character lattice into "
        "word lattices, by expanding the subpaths between \"whitespaces\" "
        "or other delimiters.\n"
        "\n"
        "Keep in mind that the expansion has an EXPONENTIAL COST with "
        "respect the maximum output degree and length of the subpaths "
        "(i.e. O(degree^length)).\n"
        "\n"
        "In practice, the exponential cost is not an issue since the "
        "subpaths are not very long. If you are facing problems with "
        "some particular lattice, you can use two pruning mechanisms "
        "to prevent the output lattices from exploding:\n"
        "\n"
        "  a) You can prune the lattices before expanding with the "
        "--beam option.\n"
        "  b) You can set a maximum length for the expanded paths "
        "--max-length option. Any path with a subpath longer than this "
        "will not be added to the output lattice.\n"
        "\n"
        "Usage: lattice-expand-subpaths [options] label-groups "
        "lat-rspecifier lat-wspecifier\n"
        " e.g.: lattice-expand-subpaths \"3 ; 4 ; 5 6\" ark:1.lat ark:1-word.lat\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    int32 max_length = std::numeric_limits<int32>::max();
    std::string symbols_table_str = "";
    bool symbols_table_text = false;

    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("beam", &beam,
                "Pruning beam (applied after lattice scaling).");
    po.Register("symbol-table", &symbols_table_str,
                "If given, all lattices will use the same symbol table which "
                "will be written to this destination file. If not provided, "
                "each lattice contains its own symbol table.");
    po.Register("symbol-table-text", &symbols_table_text,
                "If a symbol table is given, it will write the table in text "
                "mode.");
    po.Register("max-length", &max_length,
                "Maximum length of a subpath.");
    // Register TaskSequencer options.
    TaskSequencerConfig task_sequencer_config;
    task_sequencer_config.Register(&po);
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
    label_group.ParseString(po.GetArg(1));
    const std::string lattice_in_str = po.GetArg(2);
    const std::string lattice_out_str = po.GetArg(3);

    fst::SymbolTable *symbol_table =
        symbols_table_str.empty()
        ? nullptr
        : fst::ReadOrCreateSymbolTable(symbols_table_str);

    TaskSequencer<ExpandLatticeTask<CompactLattice, CompactLatticeWriter>>
        task_sequencer(task_sequencer_config);
    SequentialCompactLatticeReader lattice_reader(lattice_in_str);
    CompactLatticeWriter lattice_writer(lattice_out_str);
    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      task_sequencer.Run(
          new ExpandLatticeTask<CompactLattice, CompactLatticeWriter>(
              lattice_reader.Key(),
              lattice_reader.Value(),
              label_group,
              acoustic_scale,
              graph_scale,
              beam,
              ExpandSubpathsOptions(max_length, false),
              &lattice_writer,
              symbol_table));

      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
    }
    task_sequencer.Wait();

    if (symbol_table) {
      KALDI_VLOG(1) << "Output symbol table contains "
                    << symbol_table->NumSymbols() << " symbols.";
      fst::WriteSymbolTable(*symbol_table,
                            symbols_table_str,
                            symbols_table_text);
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return 1;
  }
}

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
#include "fstext/fstext-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/label-group.h"
#include "fstext/make-preceding-symbols-same-class.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"

namespace kaldi {

using fst::PrecedingSymbolsSameClassOptions;
using fst::MakePrecedingSymbolsSameClass;
using fst::LabelGroup;

template<typename L, typename LW>
class MakePrecedingSymbolsSameClassTask {
 public:
  typedef typename L::Arc Arc;
  typedef typename L::Arc::Label Label;

  MakePrecedingSymbolsSameClassTask(
      const std::string &key,
      const L &lat,
      const LabelGroup<Label> &label_group,
      const bool make_following_symbols_same_class,
      LW *lattice_writer)
      : key_(key),
        lat_(lat),
        label_group_(label_group),
        lattice_writer_(lattice_writer),
        make_following_symbols_same_class_(make_following_symbols_same_class)
  {}

  virtual ~MakePrecedingSymbolsSameClassTask() {
    lattice_writer_->Write(key_, lat_);
  }

  void operator()() {
    Timer timer;

    auto orig_num_states = lat_.NumStates();
    auto orig_num_arcs = NumArcs(lat_);

    PrecedingSymbolsSameClassOptions opts;
    opts.propagate_epsilon_class = true;

    if (make_following_symbols_same_class_) {
      Invert(&lat_);
    }

    std::vector<Label> state_class;
    MakePrecedingSymbolsSameClass(label_group_, &lat_, &state_class, opts);

    if (make_following_symbols_same_class_) {
      Invert(&lat_);
    }

    KALDI_LOG << "Lattice " << key_
              << " expanded #states from " << orig_num_states
              << " to " << lat_.NumStates()
              << " and #arcs from " << orig_num_arcs
              << " to " << NumArcs(lat_)
              << " in " << timer.Elapsed() << " seconds.";
  }

 protected:
  const std::string key_;
  L lat_;
  LabelGroup<Label> label_group_;
  LW *lattice_writer_;
  bool make_following_symbols_same_class_;
};

}  // namespace kaldi

int main(int argc, char **argv) {
  try {
    using namespace kaldi;
    typedef CompactLattice::Arc::Label Label;

    const char *usage =
        "";

    std::string other_groups_str = "";
    bool make_following_symbols_same_class = false;
    ParseOptions po(usage);
    po.Register("make-following-symbols-same-class",
                &make_following_symbols_same_class,
                "If true, transform the lattice so that the outgoing symbols "
                "from each are of the same class.");
    po.Register("other-groups", &other_groups_str,
                "Specific labels to group as words. Groups are separated "
                "with a semicolon, and labels within a group are separated "
                "by spaces.");
    // Register TaskSequencer options.
    TaskSequencerConfig task_sequencer_config;
    task_sequencer_config.Register(&po);
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }


    LabelGroup<Label> label_group;
    // group of word delimiters (no subpath expansion)
    if (!label_group.ParseStringSingleGroup(po.GetArg(1))) {
      KALDI_ERR << "Invalid set of non-expandable labels: \""
                << po.GetArg(1) << "\"";
    }
    // others groups to be expanded (e.g. numbers)
    if (!label_group.ParseStringMultipleGroups(other_groups_str)) {
      KALDI_ERR << "Invalid sets of additional label groups: \""
                << other_groups_str << "\"";
    }

    typedef MakePrecedingSymbolsSameClassTask<
      CompactLattice, CompactLatticeWriter> Task;
    TaskSequencer<Task> task_sequencer(task_sequencer_config);
    SequentialCompactLatticeReader lattice_reader(po.GetArg(2));
    CompactLatticeWriter lattice_writer(po.GetArg(3));
    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      task_sequencer.Run(
          new Task(
              lattice_reader.Key(),
              lattice_reader.Value(),
              label_group,
              make_following_symbols_same_class,
              &lattice_writer));
    }
    task_sequencer.Wait();
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return 1;
  }
}

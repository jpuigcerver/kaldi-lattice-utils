// fstbin/fst-compose.cc
//
// MIT License
//
// Copyright (c) 2018 Joan Puigcerver <joapuipe@gmail.com>
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
#include "fstext/fstext-lib.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/kaldi-thread.h"

namespace kaldi {

class ComposeTask {
 public:
  ComposeTask(const std::string &key, const fst::VectorFst<fst::StdArc> &fst1,
              const fst::VectorFst<fst::StdArc> &fst2,
              const fst::StdArc::Label phi_label,
              TableWriter<fst::VectorFstHolder> *fst_writer,
              size_t *n_done, size_t *n_fail)
      : key_(key), fst1_(fst1), fst2_(fst2), phi_label_(phi_label),
        fst_writer_(fst_writer), n_done_(n_done), n_fail_(n_fail) {}

  ~ComposeTask() {
    if (fst3_.Start() == fst::kNoStateId) {
      KALDI_WARN << "Empty fst for utterance " << key_ << " (wrong LM?)";
      (*n_fail_)++;
    } else {
      fst_writer_->Write(key_, fst3_);
      (*n_done_)++;
    }
  }

  void operator()() {
    fst::ArcSort(&fst1_, fst::OLabelCompare<fst::StdArc>());
    fst::ArcSort(&fst2_, fst::ILabelCompare<fst::StdArc>());
    if (phi_label_ > 0) {
      fst::PropagateFinal(phi_label_, &fst2_);
    }
    if (phi_label_ > 0) { fst::PhiCompose(fst1_, fst2_, phi_label_, &fst3_); }
    else { fst::Compose(fst1_, fst2_, &fst3_); }
  }

 private:
  const std::string key_;
  fst::VectorFst<fst::StdArc> fst1_, fst2_, fst3_;
  const fst::StdArc::Label phi_label_;
  TableWriter<fst::VectorFstHolder> *fst_writer_;
  size_t *n_done_, *n_fail_;
};

}  // namespace kaldi

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;

    const char *usage =
        "Composes FSTs.\n"
        "\n"
        "Usage: fst-compose [options] fst-rspecifier1 "
        "(fst-rspecifier2|fst-rxfilename2) fst-wspecifier\n"
        "  e.g.: fst-compose ark:1.fsts ark:2.fsts ark:composed.fsts\n"
        "  or: fst-compose ark:1.fsts G.fst ark:composed.fsts\n";

    ParseOptions po(usage);

    int32 phi_label = fst::kNoLabel; // == -1
    po.Register("phi-label", &phi_label,
                "If >0, the label on backoff arcs of the LM.");

    // Register TaskSequencer options.
    TaskSequencerConfig task_sequencer_config;
    task_sequencer_config.Register(&po);

    po.Read(argc, argv);
    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    // e.g. 0 not allowed.
    KALDI_ASSERT(phi_label > 0 || phi_label == fst::kNoLabel);

    std::string fst_rspecifier1 = po.GetArg(1),
        arg2 = po.GetArg(2),
        fst_wspecifier = po.GetArg(3);
    size_t n_done = 0, n_fail = 0;

    SequentialTableReader<fst::VectorFstHolder> fst_reader1(fst_rspecifier1);
    TableWriter<fst::VectorFstHolder> fst_writer(fst_wspecifier);
    TaskSequencer<ComposeTask> task_sequencer(task_sequencer_config);

    if (ClassifyRspecifier(arg2, NULL, NULL) == kNoRspecifier) {
      std::unique_ptr<VectorFst<StdArc>> fst2(fst::ReadFstKaldi(arg2));
      for (; !fst_reader1.Done(); fst_reader1.Next()) {
        std::string key = fst_reader1.Key();
        KALDI_VLOG(1) << "Processing fst for key " << key;
        task_sequencer.Run(
            new ComposeTask(key, fst_reader1.Value(), *fst2, phi_label,
                            &fst_writer, &n_done, &n_fail));
      }
    } else {
      RandomAccessTableReader<fst::VectorFstHolder> fst_reader2(arg2);
      for (; !fst_reader1.Done(); fst_reader1.Next()) {
        const std::string key = fst_reader1.Key();
        KALDI_VLOG(1) << "Processing fst for key " << key;
        if (!fst_reader2.HasKey(key)) {
          KALDI_WARN << "Not producing output for utterance " << key
                     << " because not present in second table.";
          n_fail++;
          continue;
        }
        task_sequencer.Run(
            new ComposeTask(key, fst_reader1.Value(), fst_reader2.Value(key),
                            phi_label, &fst_writer, &n_done, &n_fail));
      }
    }
    // Wait for all tasks.
    task_sequencer.Wait();
    KALDI_LOG << "Done " << n_done << " fsts; failed for " << n_fail;
    return (n_done != 0 ? 0 : 1);
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

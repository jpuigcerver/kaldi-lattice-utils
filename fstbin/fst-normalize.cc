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
#include "fstext/normalize_fst.h"


int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::LogArc;

    bool use_log = true;

    const char *usage =
        "Normalize costs in FSTs.\n"
        "\n"
        "By default, the normalization is done in the log-semiring, so that\n"
        "the total sum of the costs of all paths is 0.\n"
        "\n"
        "You can also normalize in the tropical-semiring, so that lowest cost\n"
        "in the FST is 0.\n"
        "\n"
        "Usage: fst-normalize [options] fsts-rspecifier fsts-wspecifier\n"
        "  e.g: fst-normalize --use-log=false ark:fst.ark arkt:fst.norm.ark\n";

    ParseOptions po(usage);
    po.Register("use-log", &use_log, "If true, normalize in the log-semiring");

    po.Read(argc, argv);
    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    const std::string fsts_rspecifier = po.GetArg(1);
    const std::string fsts_wspecifier = po.GetArg(2);

    SequentialTableReader<fst::VectorFstHolder> fst_reader(fsts_rspecifier);
    TableWriter<fst::VectorFstHolder> fst_writer(fsts_wspecifier);

    for (; !fst_reader.Done(); fst_reader.Next()) {
      const string& fst_key = fst_reader.Key();

      if (use_log) {
        // Convert StdArc to LogArc
        VectorFst<LogArc> log_fst;
        fst::ArcMap(fst_reader.Value(), &log_fst,
                    fst::WeightConvertMapper<StdArc, LogArc>());
        // Normalize LogArc
        NormalizeFst(&log_fst);
        // Convert LogArc to StdArc
        fst::ArcMap(log_fst, &fst_reader.Value(),
                    fst::WeightConvertMapper<LogArc, StdArc>());
      } else {
        // Normalize StdArc
        NormalizeFst(&fst_reader.Value());
      }

      // Always write StdArc
      fst_writer.Write(fst_key, fst_reader.Value());
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

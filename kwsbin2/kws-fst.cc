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

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::LogArc;

    bool use_inputs = false;

    const char *usage =
        "Usage: kws-fst [options] fsts1-rspecifier fsts2-rspecifier\n"
        "  e.g: kws-fst --use-inputs=true ark:fst1.ark ark:fst2.ark\n";

    ParseOptions po(usage);
    po.Register("use-inputs", &use_inputs,
                "If true, use input symbols to compose the FSTs");

    po.Read(argc, argv);
    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    const std::string fsts_rspecifier1 = po.GetArg(1);
    const std::string fsts_rspecifier2 = po.GetArg(1);

    SequentialTableReader<fst::VectorFstHolder> fst_reader1(fsts_rspecifier1);
    for (; !fst_reader1.Done(); fst_reader1.Next()) {
      const string& key1 = fst_reader1.Key();

      VectorFst<LogArc> fst1;
      fst::ArcMap(fst_reader1.Value(), &fst1,
                  fst::WeightConvertMapper<StdArc, LogArc>());
      fst_reader1.FreeCurrent();


      if (use_inputs) {
        fst::Project(&fst1, fst::PROJECT_INPUT);
        fst::ArcSort(&fst1, fst::OLabelCompare<LogArc>());
      } else {
        fst::Project(&fst1, fst::PROJECT_OUTPUT);
        fst::ArcSort(&fst1, fst::ILabelCompare<LogArc>());
      }


      SequentialTableReader<fst::VectorFstHolder> fst_reader2(fsts_rspecifier2);
      for (; !fst_reader2.Done(); fst_reader2.Next()) {
        const string& key2 = fst_reader2.Key();

        VectorFst<LogArc> fst2;
        fst::ArcMap(fst_reader2.Value(), &fst2,
                    fst::WeightConvertMapper<StdArc, LogArc>());
        fst_reader2.FreeCurrent();

        VectorFst<LogArc> comp_fst;
        if (use_inputs) {
          fst::ArcSort(&fst2, fst::ILabelCompare<LogArc>());
          fst::Compose(fst1, fst2, &comp_fst);
        } else {
          fst::ArcSort(&fst2, fst::OLabelCompare<LogArc>());
          fst::Compose(fst2, fst1, &comp_fst);
        }

        std::vector<LogArc::Weight> costs;
        fst::ShortestDistance(comp_fst, &costs, true);

        const double logprob = -costs[0].Value();
        std::cout << key1 << " " << key2 << " "
                  << std::scientific << logprob << std::endl;
      }
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

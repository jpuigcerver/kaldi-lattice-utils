// fstbin/fst-info.cc

// Copyright 2009-2011  Microsoft Corporation
//                2013  Johns Hopkins University (author: Daniel Povey)
//                2018  Joan Puigcerver

// See ../../COPYING for clarification regarding multiple authors
//
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

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/fst-info.h"

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;

    const char *usage =
        "Prints out a summary of the FSTs in a set of archives or\n"
        "the individual information about each of FST\n"
        "(--summary=false), similar to what fstinfo does.\n"
        "Usage: fst-info [options] fst-rspecifier1 [fst-rspecifier2 ...]\n"
        " e.g.: fst-info --summary=false ark:1.fsts ark:2.fsts\n";

    ParseOptions po(usage);
    bool summary = true;
    std::string include_rxfilename;
    std::string exclude_rxfilename;

    po.Register("summary", &summary,
                "If true, summarizes the information of all FSTs.");
    po.Register("include", &include_rxfilename,
                "Text file, the first field of each "
                "line being interpreted as the "
                "utterance-id whose FST will be included");
    po.Register("exclude", &exclude_rxfilename,
                "Text file, the first field of each "
                "line being interpreted as an utterance-id "
                "whose FST will be excluded");

    po.Read(argc, argv);

    if (po.NumArgs() < 1) {
      po.PrintUsage();
      exit(1);
    }

    if (include_rxfilename != "" && exclude_rxfilename != "") {
      KALDI_ERR << "should not have both --exclude and --include option!";
    }

    typedef kaldi::SequentialTableReader<fst::VectorFstHolder> SequentialVectorFstReader;
    fst::FstSummaryAcc summary_acc;
    for (int32 a = 1; a <= po.NumArgs(); ++a) {
      const std::string fsts_rspecifier = po.GetArg(a);
      if (summary) {
        fst::UpdateFstSummary<SequentialVectorFstReader>(fsts_rspecifier,
                                                         &summary_acc);
      } else {
        fst::PrintFstInfo<SequentialVectorFstReader>(fsts_rspecifier);
      }
    }
    if (summary) {
      std::cout << summary_acc;
    }

    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}
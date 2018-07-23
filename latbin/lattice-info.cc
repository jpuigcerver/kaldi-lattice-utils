// latbin/lattice-info.cc

// Copyright 2009-2011  Microsoft Corporation
//                2013  Johns Hopkins University (author: Daniel Povey)
//                2017  Joan Puigcerver

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
#include "lat/kaldi-lattice.h"
#include "fstext/fst-info.h"
#include "fstext/label-group.h"

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using fst::LabelGroup;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    typedef fst::StdArc::Label Label;

    const char *usage =
        "Prints out a summary of the lattices in a set of archives or\n"
        "the individual information about each of the lattices\n"
        "(--summary=false), similar to what fstinfo does.\n"
        "Usage: lattice-info [options] lattice-rspecifier1 [lattice-rspecifier2 ...]\n"
        " e.g.: lattice-info --summary=false ark:1.lats ark:2.lats\n";

    ParseOptions po(usage);
    bool compact = true;
    bool summary = true;
    std::string include_rxfilename;
    std::string exclude_rxfilename;
    std::string label_groups_str;

    po.Register("compact", &compact,
                "If true, work with lattices in compact form.");
    po.Register("summary", &summary,
                "If true, summarizes the information of all lattices.");
    po.Register("include", &include_rxfilename,
                "Text file, the first field of each "
                "line being interpreted as the "
                "utterance-id whose lattices will be included");
    po.Register("exclude", &exclude_rxfilename,
                "Text file, the first field of each "
                "line being interpreted as an utterance-id "
                "whose lattices will be excluded");
    po.Register("label-groups", &label_groups_str,
                "Groups of labels to form subpaths. Groups are separated "
                "with a semicolon, and labels within a group are separated "
                "by spaces.");

    po.Read(argc, argv);
    if (po.NumArgs() < 1) {
      po.PrintUsage();
      exit(1);
    }

    LabelGroup<Label> label_group;
    if (!label_group.ParseStringMultipleGroups(label_groups_str)) {
      KALDI_ERR << "Invalid sets of label groups: \""
                << label_groups_str << "\"";
    }

    if (include_rxfilename != "" && exclude_rxfilename != "") {
      KALDI_ERR << "should not have both --exclude and --include option!";
    }


    const LabelGroup<Label>* label_group_ptr =
        label_group.NumGroups() > 1 ? &label_group : nullptr;
    fst::FstSummaryAcc summary_acc;
    for (int32 a = 1; a <= po.NumArgs(); ++a) {
      const std::string lats_rspecifier = po.GetArg(a);
      if (compact) {
        if (summary) {
          fst::UpdateFstSummary<SequentialCompactLatticeReader, LabelGroup<Label>>(
              lats_rspecifier, &summary_acc, label_group_ptr);
        } else {
          fst::PrintFstInfo<SequentialCompactLatticeReader, LabelGroup<Label>>(
              lats_rspecifier, label_group_ptr);
        }
      } else {
        if (summary) {
          fst::UpdateFstSummary<SequentialLatticeReader, LabelGroup<Label>>(
              lats_rspecifier, &summary_acc, label_group_ptr);
        } else {
          fst::PrintFstInfo<SequentialLatticeReader, LabelGroup<Label>>(
              lats_rspecifier, label_group_ptr);
        }
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

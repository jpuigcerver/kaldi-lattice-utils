// kwsbin2/utils.h

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

#ifndef KALDI_KWSBIN2_UTILS_H_
#define KALDI_KWSBIN2_UTILS_H_

#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

namespace kaldi {

template <typename LabelMap, typename Container>
void AssignGroupToLabels(
    const typename LabelMap::mapped_type& group,
    const Container& labels, LabelMap* label_group) {
  for (const auto label : labels) {
    auto r = label_group->emplace(label, group);
    if (!r.second) {
      KALDI_ERR << "Each label must be assigned to one group at most. "
                << "Label " << label << " was assigned to both groups "
                << r.first->second << " and " << group << ".";
    }
  }
}

template <typename LabelSet, typename LabelMap, typename GroupSet>
void ParseSeparatorGroups(const std::string& wspaces_str,
                          const std::string& separator_groups_str,
                          LabelSet* wspace_labels,
                          LabelMap* label_group,
                          GroupSet* group_inc_count) {
  KALDI_ASSERT(wspace_labels != nullptr);
  KALDI_ASSERT(label_group != nullptr);
  KALDI_ASSERT(group_inc_count != nullptr);

  label_group->clear();
  group_inc_count->clear();

  // Epsilon is mapped to the special group 0
  (*label_group)[0] = 0;
  // All labels not explicitly assigned to any group are assigned to a
  // special group which always must increase the count of words.
  group_inc_count->insert(
      std::numeric_limits<typename GroupSet::value_type>::max());

  // Whitespace labels, which do not increment the count of words.
  {
    std::vector<typename LabelMap::key_type> tmp;
    kaldi::SplitStringToIntegers(wspaces_str, " ", true, &tmp);
    if (tmp.empty()) {
      KALDI_ERR << "At least one label must be specified as a whitespace "
                << "separator!";
    }
    AssignGroupToLabels(1, tmp, label_group);
    wspace_labels->insert(tmp.begin(), tmp.end());
  }

  // Other groups of separators, which do count as words (e.g. dot, colon, etc).
  {
    std::vector<std::string> tmp;
    kaldi::SplitStringToVector(separator_groups_str, ";", true, &tmp);
    for (size_t i = 0; i < tmp.size(); ++i) {
      std::vector<typename LabelMap::key_type> tmp2;
      kaldi::SplitStringToIntegers(tmp[i], " ", true, &tmp2);
      AssignGroupToLabels(i + 2, tmp2, label_group);
      group_inc_count->emplace(i + 2);
    }
  }
}

}  // namespace kaldi

namespace fst {

// Convert a CompactLattice to a FST where the cost of each arc is the
// total cost of the CompactLatticeArc (acoustic + graph costs), the
// input label is the word/char in the CompactLatticeArc, and the output label
// is an integer that maps to the initial and end frames segmentation.
template <typename Arc>
int32 CompactLatticeToSegmentFst(
    const kaldi::CompactLattice &clat,
    MutableFst<Arc> *ofst,
    std::vector<std::tuple<int32, int32>> *label_to_segment) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  ofst->DeleteStates();

  // Compute the times for each state in the lattice, so that we can get
  // the segmentation of each symbol.
  std::vector<int32> times;
  const int32 total_frames = kaldi::CompactLatticeStateTimes(clat, &times);

  // Add states to the output fst.
  for (StateId s = 0; s < clat.NumStates(); ++s) {
    ofst->SetFinal(ofst->AddState(), ConvertToCost(clat.Final(s)));
  }
  ofst->SetStart(clat.Start());

  map<std::tuple<int32, int32>, int32> segm_to_label;
  segm_to_label[std::make_tuple<int32, int32>(0, 0)] = 0;

  // Add arcs to the output fst.
  for (StateIterator<kaldi::CompactLattice> siter(clat); !siter.Done();
       siter.Next()) {
    const StateId s = siter.Value();
    for (ArcIterator<kaldi::CompactLattice> aiter(clat, s); !aiter.Done();
         aiter.Next()) {
      const kaldi::CompactLatticeArc& arc = aiter.Value();
      const auto segm = std::make_tuple(times[s], times[arc.nextstate]);
      const Label new_olabel =
          segm_to_label.emplace(segm, segm_to_label.size()).first->second;
      const double new_weight = ConvertToCost(arc.weight);
      ofst->AddArc(s, Arc(arc.ilabel, new_olabel, new_weight, arc.nextstate));
    }
  }

  // Create the vector mapping from label (int32) to segment.
  if (label_to_segment != nullptr) {
    label_to_segment->resize(segm_to_label.size());
    for (const auto& kv : segm_to_label) {
      (*label_to_segment)[kv.second] = kv.first;
    }
  }

  return total_frames;
}

template<typename Arc>
void SymbolToPathSegmentationFst(
    MutableFst<Arc>* mfst,
    const std::vector <std::tuple<int32, int32>> &label_to_segm) {
  typedef typename Arc::Weight Weight;
  const auto NS = mfst->NumStates();
  for (auto s0 = mfst->Start(); s0 < NS; ++s0) {
    for (MutableArcIterator<MutableFst<Arc>> ait(mfst, s0); !ait.Done(); ait.Next()) {
      Arc arc = ait.Value();
      const auto s1 = arc.nextstate;
      // path is made of a single arc, we need at least two arcs in a subpath
      // to properly represent the start/end of a segment.
      if (s0 == mfst->Start() && mfst->Final(s1) != Weight::Zero()) {
        const auto t0 = std::get<0>(label_to_segm[arc.olabel]) + 1;
        const auto t1 = std::get<1>(label_to_segm[arc.olabel]) + 1;
        // Replace output label of the existing arc, and point to a new state
        arc.olabel = t0;
        arc.nextstate = mfst->AddState();
        // Add arc from the new state to the previous destination state
        mfst->AddArc(arc.nextstate, Arc(0, t1, Weight::One(), s1));
      } else if (s0 == mfst->Start()) {
        arc.olabel = std::get<0>(label_to_segm[arc.olabel]) + 1;
      } else if (mfst->Final(s1) != Weight::Zero()) {
        arc.olabel = std::get<1>(label_to_segm[arc.olabel]) + 1;
      } else {
        arc.olabel = 0;
      }
      ait.SetValue(arc);
    }
  }
}

}  // namespace fst

#endif  // KALDI_KWSBIN2_UTILS_H_
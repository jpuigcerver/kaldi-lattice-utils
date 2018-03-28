// kwsbin2/lattice-index.cc

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

#include <algorithm>
#include <string>

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"
#include "util/text-utils.h"
#include "fstext/fstext-utils2.h"
#include <fst/script/print.h>

namespace kaldi {

template <typename LabelMap, typename GroupSet>
size_t DisambiguateStatesByWordCount(
    CompactLattice* clat,
    const LabelMap& label_group,
    const GroupSet& group_inc_count,
    std::vector<typename LabelMap::mapped_type>* state_group,
    std::vector<size_t>* state_word_count) {
  CompactLattice tmp_lat;
  std::vector<typename LabelMap::mapped_type> state_group_tmp;
  fst::DisambiguateStatesByInputLabelGroup(
      *clat, label_group, &tmp_lat, &state_group_tmp);

  return fst::DisambiguateStatesByGroupTransitionsLength(
      tmp_lat, state_group_tmp, group_inc_count, clat, state_word_count,
      state_group);
}

}  // namespace kaldi


template <typename Arc>
size_t NumArcs(const fst::Fst<Arc>& ifst) {
  using namespace fst;
  typedef fst::Fst<Arc> Fst;
  size_t num_arcs = 0;
  for (StateIterator<Fst> sit(ifst); !sit.Done(); sit.Next()) {
    num_arcs += ifst.NumArcs(sit.Value());
  }
  return num_arcs;
}

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
  group_inc_count->insert(std::numeric_limits<size_t>::max());

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

template <typename Arc, typename LabelSet>
void RemoveWhitespaceArcs(
    fst::MutableFst<Arc>* ifst,
    const std::vector< std::tuple<typename Arc::Label, size_t> >& ilabel_info,
    const LabelSet& wspace_labels) {
  KALDI_ASSERT(ifst != nullptr);
  typedef fst::MutableFst<Arc> FST;

  auto dummyFinal = ifst->AddState();
  for (fst::StateIterator<FST> sit(*ifst); !sit.Done(); sit.Next()) {
    const auto u = sit.Value();
    for (fst::MutableArcIterator<FST> ait(ifst, u); !ait.Done(); ait.Next()) {
      auto arc = ait.Value();
      KALDI_ASSERT(arc.ilabel <= ilabel_info.size());
      const auto label = std::get<0>(ilabel_info[arc.ilabel]);
      if (wspace_labels.count(label)) {
        arc.nextstate = dummyFinal;
        ait.SetValue(arc);
      }
    }
  }
  fst::Connect(ifst);
}

namespace kaldi {

template <typename Arc, typename Int1, typename Int2>
void CompactLatticeToWordCountSegmFst(
    const CompactLattice& clat,
    const std::vector<Int1>& state_times,
    const std::vector<Int2>& state_word_count,
    fst::MutableFst<Arc>* ofst,
    std::vector< std::tuple<typename Arc::Label, Int2> >* ilabels,
    std::vector< std::tuple<Int1, Int1> >* olabels) {
  using namespace fst;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  ofst->DeleteStates();
  ilabels->clear();
  olabels->clear();
  if (clat.Start() == kNoStateId) return;

  KALDI_ASSERT(state_times.size() == clat.NumStates());
  KALDI_ASSERT(state_word_count.size() == clat.NumStates());

  // Add states to the output fst and set final weight.
  for (StateId s = 0; s < clat.NumStates(); ++s) {
    ofst->SetFinal(ofst->AddState(),
                   clat.Final(s).Weight().Value1() +
                   clat.Final(s).Weight().Value2());
  }
  // Set the initial state for the
  ofst->SetStart(clat.Start());

  // Use these to map from tuples to label IDs.
  std::map<std::tuple<Label, Int2>, Label> ilabels_map;
  std::map<std::tuple<Int1, Int1>, Label> olabels_map;
  ilabels_map[std::make_tuple(0, 0)] = 0;
  olabels_map[std::make_tuple(0, 0)] = 0;

  for (StateIterator<CompactLattice> sit(clat); !sit.Done(); sit.Next()) {
    const StateId s = sit.Value();
    for (ArcIterator<CompactLattice> ait(clat, s); !ait.Done(); ait.Next()) {
      const auto& arc = ait.Value();
      const auto new_weight =
          arc.weight.Weight().Value1() + arc.weight.Weight().Value2();
      // Get input label for the equivalent arc in the output fst.
      const auto ituple =
          std::make_tuple(arc.ilabel, state_word_count[arc.nextstate]);
      const auto ilabel =
          ilabels_map.emplace(ituple, ilabels_map.size()).first->second;
      // Get output label for the equivalent arc in the output fst.
      const auto otuple =
          std::make_tuple(state_times[s], state_times[arc.nextstate]);
      const auto olabel =
          olabels_map.emplace(otuple, olabels_map.size()).first->second;
      // Add arc to the output fst.
      ofst->AddArc(s, Arc(ilabel, olabel, new_weight, arc.nextstate));
    }
  }

  ilabels->resize(ilabels_map.size());
  for (const auto& kv : ilabels_map) {
    (*ilabels)[kv.second] = kv.first;
  }

  olabels->resize(olabels_map.size());
  for (const auto& kv : olabels_map) {
    (*olabels)[kv.second] = kv.first;
  }

  auto outprops =
      clat.Properties(kFstProperties & ~(kAcceptor | kNotAcceptor), false);
  outprops |= kNotAcceptor;
  ofst->SetProperties(outprops, kFstProperties);
}


}  // namespace kaldi







int main(int argc, char** argv) {
  try {
    using namespace kaldi;
    using namespace fst;
    typedef CompactLatticeArc::Label Label;

    const char* usage =
        "Usage: \n"
        "  lattice-positional-index-from-chars <wspace-labels> "
        "<separator-label-groups> <lattice-ark>";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat delta = fst::kDelta;
    int nbest = 100;

     po.Register("delta", &delta, "Tolerance used in determinization. "
                "The smaller the better.");
    po.Register("acoustic-scale", &acoustic_scale,
                 "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Register("nbest", &nbest, "Extract this number of n-best hypothesis.");
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    // Parse separator symbols from arguments
    std::unordered_set<Label> wspace_labels;
    std::unordered_map<Label, size_t> label_group;
    std::unordered_set<size_t> groups_inc_word_count;
    ParseSeparatorGroups(po.GetArg(1), po.GetArg(2), &wspace_labels,
                         &label_group, &groups_inc_word_count);

    // Scaling scores
    std::vector<std::vector<double> > scale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = graph_scale;
    scale[1][1] = acoustic_scale;

    const std::string lattice_rspecifier = po.GetArg(3);
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(scale, &clat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(insertion_penalty, &clat);
      // Lattice prunning
      if (beam != std::numeric_limits<BaseFloat>::infinity())
        PruneLattice(beam, &clat);
      // If needed, sort the compact lattice in topological order
      TopSortCompactLatticeIfNeeded(&clat);
      // Remove epsilon arcs from the compact lattice.
      RmEpsilon(&clat);
      // Make sure that all sequences arriving to each state have the same
      // WORD count (even if the number of characters is different).
      std::vector<size_t> state_group;
      std::vector<size_t> state_word_count;
      const size_t max_word_count = DisambiguateStatesByWordCount(
          &clat, label_group, groups_inc_word_count,
          &state_group, &state_word_count);

      // Compute forward and backward scores for each state.
      std::vector<double> fw, bw;
      const double total_cost =
          -ComputeLatticeAlphasAndBetas(clat, false, &fw, &bw);
      for (auto& v : fw) { v = -v; }
      for (auto& v : bw) { v = -v; }

      // Compute times for each state.
      std::vector<int32> state_times;
      const auto num_frames = kaldi::CompactLatticeStateTimes(
          clat, &state_times);

      VectorFst<LogArc> log_fst;
      std::vector< std::tuple<Label, size_t> > ilabels;
      std::vector< std::tuple<int32, int32> > olabels;
      CompactLatticeToWordCountSegmFst(clat, state_times, state_word_count,
                                       &log_fst, &ilabels, &olabels);
      GroupFactorFst(&log_fst, fw, bw, state_group);

      // Remove paths from the fst which represent groups of whitespaces,
      // since we don't want to index whitespaces.
      RemoveWhitespaceArcs(&log_fst, ilabels, wspace_labels);

      using AM1 = WeightConvertMapper<LogArc, StdArc>;
      using AM2 = RmWeightMapper<StdArc, StdArc>;

      // Determinize in the log semiring to sum weights of all paths of the
      // same (word, position).
      DeterminizeFstOptions<LogArc> dopts1(delta, 0, DETERMINIZE_FUNCTIONAL);
      ArcMapFst<LogArc, StdArc, AM1> det_fst1(
          DeterminizeFst<LogArc>(
              ProjectFst<LogArc>(log_fst, PROJECT_INPUT),
              dopts1),
          AM1());

      // Determinize in the tropical semiring to obtain the best segmentation
      // for each (word, position).
      DeterminizeFstOptions<StdArc> dopts2(delta, 0, DETERMINIZE_DISAMBIGUATE);
      ArcMapFst<StdArc, StdArc, AM2> det_fst2(
          DeterminizeFst<StdArc>(
              ArcMapFst<LogArc, StdArc, AM1>(log_fst, AM1()),
              dopts2),
          AM2());

      VectorFst<StdArc> nbest_fst;
      ShortestPath(ComposeFst<StdArc>(det_fst1, det_fst2), &nbest_fst, nbest);


      std::vector< VectorFst<StdArc> > nbest_fsts;
      ConvertNbestToVector(nbest_fst, &nbest_fsts);
      for (const auto& nbfst : nbest_fsts) {
        std::vector<StdArc::Label> nb_isyms, nb_osyms;
        StdArc::Weight nb_cost;
        GetLinearSymbolSequence(nbfst, &nb_isyms, &nb_osyms, &nb_cost);
        if (nb_isyms.empty() || nb_osyms.empty()) {
          KALDI_WARN << "Ignoring empty word in lattice \""
                     << lattice_key << "\"";
          continue;
        }

        const auto word_pos = std::get<1>(ilabels[nb_isyms.front()]);
        std::cout << lattice_key << " " << total_cost - nb_cost.Value() << " " << word_pos;
        std::cout << " " << std::get<0>(olabels[nb_osyms.front()])
                  << " " << std::get<1>(olabels[nb_osyms.back()]);
        for (auto s : nb_isyms) {
          KALDI_ASSERT(std::get<1>(ilabels[s]) == word_pos);
          std::cout << " " << std::get<0>(ilabels[s]);
        }
        std::cout << std::endl;
      }




      /*
      std::cout << lattice_key << " ; max_len = " << max_len << " ; "
                << "NumStates = " << clat.NumStates() << " vs. "
                << clat_tmp.NumStates() << " ; "
                << "NumArcs = " << NumArcs(clat) << " vs. "
                << NumArcs(clat_tmp) << std::endl;
      */
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
}

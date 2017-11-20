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

namespace fst {

template <typename Arc>
size_t DisambiguateStateInputSequenceLengthFromCharacters(
    const Fst<Arc>& ifst,
    const std::vector<std::vector<typename Arc::Label>>& separator_groups,
    const std::vector<bool>& separator_group_inc_length,
    MutableFst<Arc>* ofst,
    std::vector<size_t>* state_input_length = nullptr,
    bool use_input = false) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(ofst != nullptr);
  KALDI_ASSERT(separator_groups.size() == separator_group_inc_length.size());
  if (ifst.Properties(kAcyclic, true) != kAcyclic) {
    KALDI_ERR << "Only acyclic transducers are supported.";
  }
  // Create map of separators.
  std::unordered_map<Label, size_t> separators;
  separators[0] = 0;
  for (size_t i = 0; i < separator_groups.size(); ++i) {
    for (const Label& a : separator_groups[i]) {
      if (!separators.emplace(a, i + 1).second) {
        KALDI_ERR << "Label " << a << " was assigned to two separator groups.";
      }
    }
  }
  if (separators.size() == 0) {
    KALDI_ERR << "At least one separator label is required.";
  }

  // Useful lambda function to obtain the group of a given label.
  auto get_label_group = [&](const Label& lbl) -> size_t {
    auto it = separators.find(lbl);
    return it != separators.end() ? it->second : separator_groups.size() + 1;
  };
  auto does_group_inc_length = [&](const size_t& group) -> bool {
    if (group == 0) return false;
    if (group > separator_group_inc_length.size()) return true;
    return separator_group_inc_length[group - 1];
  };

  ofst->DeleteStates();
  if (state_input_length != nullptr) state_input_length->clear();
  if (ifst.Start() < 0) return 0;

  // In the new fst, the states are represented by tuples:
  // (input length, state, group) from the old fst.
  // This map is used to map these tuples to StateId in the output fst.
  std::map<std::tuple<size_t, StateId, size_t>, StateId> state_map;
  // This queue is used to build the output fst.
  std::queue<std::tuple<size_t, StateId, size_t>> Q;

  // We do a first pass through the fst in order to get all the states of the
  // new fst. We need to do this first pass in order to know all the states
  // before creating any arc, and thus, creating the new fst to be sorted in
  // topological order.
  state_map[std::make_tuple(0, ifst.Start(), 0)] = -1;
  Q.push(std::make_tuple(0, ifst.Start(), 0));
  size_t max_len = 0;
  while (!Q.empty()) {
    const size_t len = std::get<0>(Q.front());
    const StateId u  = std::get<1>(Q.front());
    const size_t pg  = std::get<2>(Q.front());
    Q.pop();
    // Update the maximum length seen so far.
    if (max_len < len) max_len = len;
    for (ArcIterator<Fst<Arc>> aiter(ifst, u); !aiter.Done(); aiter.Next()) {
      const auto arc = aiter.Value();
      const auto lbl = use_input ? arc.ilabel : arc.olabel;
      const auto g = get_label_group(lbl);
      // Only when the label group should increase the length and we have
      // changed the group from the previous node, we increase the length.
      // Otherwise, keep the previous length.
      const size_t next_len =
          (does_group_inc_length(g) && g != pg) ? (len + 1) : (len);
      // Create a tuple representing a new state. The group of the new state
      // is the group of label, if it is non-epsilon. If the label is epsilon,
      // it propagates the label from the previous state.
      const auto t = std::make_tuple(next_len, arc.nextstate, lbl ? g : pg);
      if (state_map.emplace(t, -1).second) {
        Q.push(t);
      }
    }
  }

  // Create states in the new fst.
  // Note: We set the NEGATIVE StateId of the new fst into state_map.
  // We will change this to a positive value when the arcs from this new state
  // have been processed.
  for (auto it = state_map.begin(); it != state_map.end(); ++it) {
    it->second = -ofst->AddState();
  }
  ofst->SetStart(0);

  // Set the state_input_length vector, if given.
  if (state_input_length != nullptr) {
    state_input_length->reserve(state_map.size());
    for (auto it = state_map.begin(); it != state_map.end(); ++it) {
      state_input_length->push_back(std::get<0>(it->first));
    }
  }

  // Traverse the fst again in order to add the arcs in the new fst.
  Q.push(std::make_tuple(0, 0, ifst.Start()));
  while (!Q.empty()) {
    const size_t len = std::get<0>(Q.front());  // Input length to new state
    const StateId u  = std::get<1>(Q.front());  // StateId in the old fst
    const Label pg   = std::get<2>(Q.front());  // Group of the state
    const StateId u2 = state_map[Q.front()];    // StateId in the new fst
    Q.pop();

    // Set final weight for the new state in the output fst.
    ofst->SetFinal(u2, ifst.Final(u));

    // Traverse edges from state u
    for (ArcIterator< Fst<Arc> > aiter(ifst, u); !aiter.Done(); aiter.Next()) {
      auto arc = aiter.Value();
      const auto lbl = use_input ? arc.ilabel : arc.olabel;
      const auto g = get_label_group(lbl);
      const size_t next_len =
          (does_group_inc_length(g) && g != pg) ? (len + 1) : (len);
      const auto t = std::make_tuple(next_len, arc.nextstate, lbl ? g : pg);
      auto state_it = state_map.find(t);  // Note: we are sure that this exists.
      if (state_it->second < 0 ) {
        state_it->second = -state_it->second;
        Q.push(t);
      }
      // Add arc to the output fst
      arc.nextstate = state_it->second;
      ofst->AddArc(u2, arc);
      if (arc.nextstate <= u2) {
        KALDI_ERR
            << "Added arc ("
            << len << ", " << u << ", " << pg << ") [" << u2 << "] -> ("
            << std::get<0>(t) << ", " << std::get<1>(t) << ", "
            << std::get<2>(t) << ") [" << arc.nextstate << "]";
      }
    }
  }

  return max_len;
}

// Given a fst that represents sequences of characters in its paths, creates
// a fst that accepts words (determined by sequences of characters in between
// special labels, separators) which where part of the original fst.
template <
  typename Arc,
  typename LabelMap = std::unordered_map<typename Arc::Label, size_t> >
void WordsFst(fst::MutableFst<Arc>* fst,
              const LabelMap& label_to_group
              bool use_input = false) {
  typedef fst::MutableFst<Arc> Fst;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(fst != NULL);
  if (fst->Start() == fst::kNoStateId) return;

  // Compute forward and backward scores for each state
  std::vector<Weight> fw, bw;
  ShortestDistance<Arc>(*fst, &fw, false);
  ShortestDistance<Arc>(*fst, &bw, true);
  const float total_cost = bw[fst->Start()].Value();

  std::vector<typename LabelMap::value_type> state_group(fw.size(), 0);
  for (StateIterator<Fst> siter(*fst); !siter.Done(); siter.Next()) {
    for (ArcIterator<Fst> aiter(*fst, siter.Value()); !aiter.Done();
         aiter.Next()) {
    }
  }

  // New final state
  const StateId sFinal = fst->AddState();

  // Convert fst to accept words in the original fst
  std::vector<Arc> new_arcs_from_init;
  for (fst::StateIterator<Fst> siter(*fst); !siter.Done(); siter.Next()) {
    const StateId s = siter.Value();
    if (s == sFinal) continue;
    std::vector<Arc> new_arcs;  // New arcs from state s
    // Add arc to the new (unique) final state, and make s not final
    if (fst->Final(s) != Weight::Zero()) {
      new_arcs.push_back(Arc(0, 0, fst->Final(s), sFinal));
      fst->SetFinal(s, Weight::Zero());
    }
    // Traverse current arcs and remove arcs with separator labels.
    // For each remove arc two epsilon arcs are added:
    //   - one from the current state (s) to the final state with cost =
    //     arc.weight * bw[arc.nextstate]. This is because the node s is
    //     the final node for a word.
    //   - one from the initial state to arc.nextstate, with cost =
    //     arc.weight * forward[s]. This is because arc.nextstate is the
    //     start of a new word.
    for (fst::ArcIterator<Fst> aiter(*fst, s); !aiter.Done(); aiter.Next()) {
      const Arc& arc = aiter.Value();
      if (separators.count(arc.ilabel)) {
        new_arcs.push_back(
            Arc(0, 0, fst::Times(arc.weight, bw[arc.nextstate]), sFinal));
        new_arcs_from_init.push_back(
            Arc(0, 0, fst::Times(arc.weight, fw[s]), arc.nextstate));

      } else {
        new_arcs.push_back(arc);
      }
    }
    // Delete all arcs from state s
    fst->DeleteArcs(s);
    // Add new arcs from state s
    for (const Arc& arc : new_arcs) {
      fst->AddArc(s, arc);
    }
  }
  // Add missing arcs from the initial state
  for (const Arc& arc: new_arcs_from_init) {
    fst->AddArc(fst->Start(), arc);
  }
  // Final cost = -total_cost, so that paths are normalized in -log [0, 1]
  fst->SetFinal(sFinal, Weight(-total_cost));
  // Remove epsilon symbols O(V^2 + V * E)
  RmEpsilon(fst);
  // Remove unnecessary states/arcs
  Connect(fst);
  // Empty strings are not allowed
  if (fst->Final(fst->Start()) != Weight::Zero()) {
    fst->SetFinal(fst->Start(), Weight::Zero());
  }
  // Push weights to the toward the initial state. This speeds up n-best list
  // retrieval.
  Push<Arc>(fst, fst::REWEIGHT_TO_INITIAL, fst::kPushWeights);
}

}  // namespace fst

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

template <typename I>
void ParseSeparatorGroups(const std::string& str,
                          std::vector<std::vector<I>>* groups,
                          std::vector<bool>* groups_inc_length) {
  groups->clear();
  groups_inc_length->clear();

  std::vector<std::string> groups_str;
  kaldi::SplitStringToVector(str, ";", true, &groups_str);
  for (const std::string& g_str : groups_str) {
    groups->push_back(std::vector<I>());
    kaldi::SplitStringToIntegers(g_str, " ", true, &groups->back());
    if (groups->back().size() < 2) {
      KALDI_ERR << "";
    }
    groups_inc_length->push_back(static_cast<bool>(groups->back().back()));
    groups->back().pop_back();
  }
}

int main(int argc, char** argv) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage = "";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    po.Register("acoustic-scale", &acoustic_scale,
                 "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Read(argc, argv);

    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    // Parse separator symbols from arguments
    std::vector< std::vector<int32> > separator_groups;
    std::vector<bool> separator_groups_inc_length;
    ParseSeparatorGroups(
        po.GetArg(1), &separator_groups, &separator_groups_inc_length);

    // Scaling scores
    std::vector<std::vector<double> > scale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = graph_scale;
    scale[1][1] = acoustic_scale;

    const std::string lattice_rspecifier = po.GetArg(2);
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
      // Make sure that all sequences arriving to each state have the same
      // WORD length (even if the number of characters is different).
      CompactLattice clat_tmp;
      std::vector<size_t> state_input_length;
      const auto max_len = DisambiguateStateInputSequenceLengthFromCharacters(
          clat, separator_groups, separator_groups_inc_length,
          &clat_tmp, &state_input_length, false);



      std::cout << lattice_key << " ; max_len = " << max_len << " ; "
                << "NumStates = " << clat.NumStates() << " vs. "
                << clat_tmp.NumStates() << " ; "
                << "NumArcs = " << NumArcs(clat) << " vs. "
                << NumArcs(clat_tmp) << std::endl;
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
}

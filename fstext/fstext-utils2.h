// fstext/fstext-utils2.h

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

#ifndef KALDI_FSTEXT_FSTEXT_UTILS2_H_
#define KALDI_FSTEXT_FSTEXT_UTILS2_H_

#include <map>
#include <queue>
#include <vector>

#include <fst/fstlib.h>

#include "base/kaldi-common.h"

namespace fst {

// Given an acyclic fst, creates an equivalent fst such that all sequences
// arriving to any state have the same length.
// The function also returns the maximum sequence length in the transducer.
//
// The argument use_input can be used to disambiguate using the input/output
// sequence length (by default, it uses the output sequences).
//
// If given, the vector state_input_length will contain for each output state
// the associated length of the sequences entering this state.
//
// The cost of this algorithm is O((V + E) * I_L), where I_L is the maximum
// number of different lengths arriving to a state. The spatial cost is the
// same.
template <typename Arc>
size_t DisambiguateStateInputSequenceLength(
    const Fst<Arc>& ifst, MutableFst<Arc>* ofst,
    std::vector<size_t>* state_input_length = nullptr,
    bool use_input = false) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(ofst != nullptr);
  if (ifst.Properties(kAcyclic, true) != kAcyclic) {
    KALDI_ERR << "DisambiguateStateInputSequenceLength() only supports acyclic "
              << "transducers.";
  }

  ofst->DeleteStates();
  if (state_input_length != nullptr) state_input_length->clear();
  if (ifst.Start() < 0) return 0;

  // In the new fst, the states are represented by tuples:
  // (input length, state) from the old fst.
  // This map is used to map these tuples to StateId in the output fst.
  std::map<std::tuple<size_t, StateId>, StateId> state_map;
  // This queue is used to build the output fst.
  std::queue<std::tuple<size_t, StateId>> Q;

  // We do a first pass through the fst in order to get all the states of the
  // new fst. We need to do this first pass in order to know all the states
  // before creating any arc, and thus, creating the new fst to be sorted in
  // topological order.
  state_map[std::make_tuple(0, ifst.Start())] = -1;
  Q.push(std::make_tuple(0, ifst.Start()));
  size_t max_len = 0;
  while (!Q.empty()) {
    const size_t len = std::get<0>(Q.front());
    const StateId u  = std::get<1>(Q.front());
    Q.pop();
    // Update the maximum length seen so far.
    if (max_len < len) max_len = len;
    for (ArcIterator<Fst<Arc>> aiter(ifst, u); !aiter.Done(); aiter.Next()) {
      auto arc = aiter.Value();
      // If the arc is labeled with an epsilon, then the length of the
      // sequence is no incresead. Otherwise, it is increased.
      const size_t next_len =
          ((use_input ? arc.ilabel : arc.olabel) == 0) ? (len) : (len + 1);
      const auto t = std::make_tuple(next_len, arc.nextstate);
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
  Q.push(std::make_tuple(0, ifst.Start()));
  while (!Q.empty()) {
    const size_t len = std::get<0>(Q.front());  // Input length to new state
    const StateId u  = std::get<1>(Q.front());  // StateId in the old fst
    const StateId u2 = state_map[Q.front()];    // StateId in the new fst
    Q.pop();

    // Set final weight for the new state in the output fst.
    ofst->SetFinal(u2, ifst.Final(u));

    // Traverse edges from state u
    for (ArcIterator< Fst<Arc> > aiter(ifst, u); !aiter.Done(); aiter.Next()) {
      auto arc = aiter.Value();
      // If the arc is labeled with an epsilon, then the length of the sequence
      // is no incresead. Otherwise, it is increased.
      const size_t next_len =
          ((use_input ? arc.ilabel : arc.olabel) == 0) ? (len) : (len + 1);
      const auto state_tuple = std::make_tuple(next_len, arc.nextstate);
      auto state_it = state_map.find(state_tuple);
      if (state_it->second < 0 ) {
        state_it->second = -state_it->second;
        Q.push(state_tuple);
      }
      arc.nextstate = state_it->second;
      // Add arc to the output fst
      ofst->AddArc(u2, arc);

      if (arc.nextstate <= u2) {
        KALDI_ERR
            << "Added arc (" << len << ", " << u << ") [" << u2 << "] -> ("
            << std::get<0>(state_tuple) << ", " << std::get<1>(state_tuple)
            << ") [" << arc.nextstate << "]";
      }
    }
  }

  return max_len;
}

// O(V)
template <typename Arc>
void AddSequenceLengthDismabiguationSymbol(
    MutableFst<Arc>* fst, std::vector<size_t>* state_input_length,
    typename Arc::Label dis_label = kNoLabel) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(fst != nullptr);
  KALDI_ASSERT(state_input_length != nullptr);
  if (fst->NumStates() != state_input_length->size()) {
    KALDI_ERR << "AddSequenceLengthDisambiguationSymbol() expects a vector "
              << "with the length of the input sequences to each state of the "
              << "FST. The length of the vector is different to the number of "
              << "states (" << state_input_length->size() << " vs. "
              << fst->NumStates() << ")";
  }
  if (fst->NumStates() == 0) return;

  // Get the maximum length amongst all sequences in the FST.
  const size_t max_length = *std::max_element(state_input_length->begin(),
                                              state_input_length->end());
  // New vector of states. State k represents that k symbols are missing so
  // that the sequences entering this state have exactly max_length symbols.
  std::vector<StateId> aux_states(max_length + 1);
  for (size_t k = 0; k <= max_length; ++k) {
    aux_states[k] = fst->AddState();
  }
  fst->SetFinal(aux_states[max_length], Weight::One());

  // We just added new states, so we need to resize and update the
  // state_input_length vector, with the input length of each of them.
  // We also need to add edges for this states, so that an arc goes from
  // state k to state k - 1 with the special disambiguation label and weight
  // equal to 1.
  state_input_length->reserve(state_input_length->size() + max_length + 1);
  for (size_t k = 0; k <= max_length; ++k) {
    state_input_length->push_back(k);
    if (k < max_length) {
      const Arc new_arc(dis_label, dis_label, Weight::One(), aux_states[k + 1]);
      fst->AddArc(aux_states[k], new_arc);
    }
  }

  // Replace final states with arcs to one of the aux_states.
  for (StateIterator< Fst<Arc> > siter(*fst); !siter.Done(); siter.Next()) {
    const StateId u = siter.Value();
    if (u >= aux_states[0]) continue;  // Skip auxiliar states.
    const Weight final_w = fst->Final(u);
    if (final_w != Weight::Zero()) {
      fst->SetFinal(u, Weight::Zero());
      fst->AddArc(u, Arc(0, 0, final_w, aux_states[(*state_input_length)[u]]));
    }
  }
}

// Given a input FST, obtain an equivalent fst with the additional property
// that the label of all input arcs to each state are part of the same group.
// The group of a label is defined by the label map label_group.
// Optionally, it can return the group of each state in the output fst and
// you can specify wether to consider input or output labels from the input fst.
template <
  typename Arc,
  typename LabelMap = std::unordered_map<typename Arc::Label, size_t> >
void DisambiguateStatesByInputLabelGroup(
    const Fst<Arc>& ifst, const LabelMap& label_group, MutableFst<Arc>* ofst,
    std::vector<typename LabelMap::mapped_type>* state_group = nullptr,
    bool use_input = false) {
  typedef Fst<Arc> FST;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Label Label;

  ofst->DeleteStates();
  if (state_group) state_group->clear();
  if (ifst.Start() == kNoStateId) {
    return;
  }

  // This map is used to map from tuples of (input_state_id, group) to
  // output_state_id.
  std::map<std::tuple<StateId, Label>, StateId> state_map;
  state_map[std::make_tuple(ifst.Start(), 0)] = 0;

  // First compute the number of states for the output fst.
  for (StateIterator<FST> sit(ifst); !sit.Done(); sit.Next()) {
    for (ArcIterator<FST> ait(ifst, sit.Value()); !ait.Done(); ait.Next()) {
      const auto arc = ait.Value();
      const auto label = use_input ? arc.ilabel : arc.olabel;
      const auto git = label_group.find(label);
      const auto grp = git != label_group.end()
          ? git->second
          : std::numeric_limits<typename LabelMap::mapped_type>::max();
      state_map.emplace(std::make_tuple(arc.nextstate, grp), -1);
    }
  }

  // Add states in the output fst.
  if (state_group) { state_group->resize(state_map.size(), 0); }
  for (auto it = state_map.begin(); it != state_map.end(); ++it) {
    it->second = ofst->AddState();
    // Write the state group in the state_group vector, if given.
    if (state_group) { (*state_group)[it->second] = std::get<1>(it->first); }
  }
  ofst->SetStart(0);

  // Add arcs in the output fst.
  for (auto it = state_map.begin(); it != state_map.end(); ++it) {
    const StateId s1 = std::get<0>(it->first);  // Source ID in the input fst
    const StateId s2 = it->second;              // Source ID in the output fst
    // Set final weight for the node in the output fst.
    ofst->SetFinal(s2, ifst.Final(s1));
    // Add arcs to the node in the output fst.
    for (ArcIterator<FST> ait(ifst, s1); !ait.Done(); ait.Next()) {
      auto arc = ait.Value();
      const auto label = use_input ? arc.ilabel : arc.olabel;
      const auto git = label_group.find(label);
      const auto grp = git != label_group.end()
          ? git->second
          : std::numeric_limits<typename LabelMap::mapped_type>::max();
      // Replace the target ID with the output fst ID
      arc.nextstate =
          state_map.find(std::make_tuple(arc.nextstate, grp))->second;
      ofst->AddArc(s2, arc);
    }
  }

  auto outprops = ifst.Properties(kFstProperties, false);
  ofst->SetProperties(outprops, kFstProperties);
}

// Given a input fst such that the label of all input arcs to each state are
// part of the same group (defined by the map label_group), this function
// will return a vector with the group of each state.
// If the input fst has states with input arcs of multiple groups, the
// function will return false, to indicate that the input fst is not correct.
template <
  typename Arc,
  typename LabelMap = std::unordered_map<typename Arc::Label, size_t> >
bool GetStatesInputLabelGroup(
    const Fst<Arc>& ifst, const LabelMap& label_group,
    std::vector<typename LabelMap::mapped_type>* state_group,
    bool use_input = false) {
  typedef Fst<Arc> FST;

  // Initialize the group of all states.
  std::vector<bool> fixed_state;
  state_group->clear();
  for (StateIterator<FST> sit(ifst); !sit.Done(); sit.Next()) {
    state_group->push_back(0);
    fixed_state.push_back(false);
  }

  for (StateIterator<FST> sit(ifst); !sit.Done(); sit.Next()) {
    for (ArcIterator<FST> ait(ifst, sit.Value()); !ait.Done(); ait.Next()) {
      const auto arc = ait.Value();
      const auto label = use_input ? arc.ilabel : arc.olabel;
      // Get the label's group.
      const auto it = label_group.find(label);
      const auto gr = it != label_group.end()
          ? it->second
          : std::numeric_limits<typename LabelMap::mapped_type>::max();
      // The state's group is the same as the grup of its input labels.
      // Thus, we need that all labels that arrive to a node are part of the
      // same group. If this is not the case, the function terminates with an
      // error.
      if (!fixed_state[arc.nextstate]) {
        fixed_state[arc.nextstate] = true;
        (*state_group)[arc.nextstate] = gr;
      } else if ((*state_group)[arc.nextstate] != gr) {
        return false;
      }
    }
  }

  return true;
}

template <typename Arc, typename GroupSet, typename Int>
Int DisambiguateStatesByGroupTransitionsLength(
    const Fst<Arc>& ifst, const std::vector<Int>& input_state_group,
    const GroupSet& group_inc_length, MutableFst<Arc>* ofst,
    std::vector<Int>* state_input_num_transition = nullptr,
    std::vector<Int>* output_state_group = nullptr) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef std::tuple<Int, StateId> StateTuple;
  if (ifst.Properties(kAcyclic, true) != kAcyclic) {
    KALDI_ERR << "Only acyclic transducers are supported.";
  }

  ofst->DeleteStates();
  if (state_input_num_transition) { state_input_num_transition->clear(); }
  if (output_state_group) { output_state_group->clear(); }
  if (ifst.Start() < 0) return 0;

  std::map<StateTuple, StateId> state_map;
  std::queue<StateTuple> Q;

  state_map[std::make_tuple(0, ifst.Start())] = -1;
  Q.push(std::make_tuple(0, ifst.Start()));
  Int max_num_transitions = 0;

  while (!Q.empty()) {
    const auto n = std::get<0>(Q.front());
    const auto u = std::get<1>(Q.front());
    KALDI_ASSERT(u < input_state_group.size());
    const auto ug = input_state_group[u];
    Q.pop();

    // Update the maximum length seen so far.
    if (max_num_transitions < n) max_num_transitions = n;

    for (ArcIterator<Fst<Arc>> ait(ifst, u); !ait.Done(); ait.Next()) {
      const auto v = ait.Value().nextstate;
      KALDI_ASSERT(v < input_state_group.size());
      const auto vg = input_state_group[v];
      const auto vn = (ug != vg && group_inc_length.count(vg)) ? (n + 1) : (n);
      const auto t = std::make_tuple(vn, v);
      if (state_map.emplace(t, -1).second) {
        Q.push(t);
      }
    }
  }

  // Reserve size to store the number of group transitions in the paths that
  // arrive to each of the output FST states.
  if (state_input_num_transition) {
    state_input_num_transition->reserve(state_map.size());
  }

  if (output_state_group) {
    output_state_group->reserve(state_map.size());
  }

  // Create states in the new fst.
  // Note: We set the NEGATIVE StateId of the new fst into state_map.
  // We will change this to a positive value when the arcs from this new state
  // have been processed.
  for (auto it = state_map.begin(); it != state_map.end(); ++it) {
    it->second = ofst->AddState();
    if (state_input_num_transition) {
      state_input_num_transition->push_back(std::get<0>(it->first));
    }
    if (output_state_group) {
      output_state_group->push_back(input_state_group[std::get<1>(it->first)]);
    }
  }
  ofst->SetStart(0);


  // Add arcs in the output fst.
  for (auto it = state_map.begin(); it != state_map.end(); ++it) {
    const auto n  = std::get<0>(it->first);
    const auto u1 = std::get<1>(it->first);
    const auto u2 = it->second;              // Source ID in the output fst
    KALDI_ASSERT(u1 < input_state_group.size());
    const auto ug = input_state_group[u1];
    // Set final weight for the node in the output fst.
    ofst->SetFinal(u2, ifst.Final(u1));
    // Add arcs to the node in the output fst.
    for (ArcIterator<Fst<Arc>> ait(ifst, u1); !ait.Done(); ait.Next()) {
      auto arc = ait.Value();
      const auto v = arc.nextstate;
      KALDI_ASSERT(v < input_state_group.size());
      const auto vg = input_state_group[v];
      const auto vn = (ug != vg && group_inc_length.count(vg)) ? (n + 1) : (n);
      const auto t = std::make_tuple(vn, v);
      arc.nextstate = state_map.find(t)->second;
      ofst->AddArc(u2, arc);
    }
  }

  auto outprops = ifst.Properties(kFstProperties, false);
  ofst->SetProperties(outprops, kFstProperties);

  return max_num_transitions;
}

// Given a state-grouped FST (a regular FST whose states are part of a given
// group), and the Forward and Backward costs of each state, transforms the
// given FST, so that its full-paths contain all subpaths from the original
// FST that traverse states of the same group.
template <typename Arc, typename Int, typename W>
void GroupFactorFst(
    MutableFst<Arc>* fst, const std::vector<Int>& state_group,
    const std::vector<W>& fw, const std::vector<W>& bw) {
  typedef MutableFst<Arc> Fst;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(fst != nullptr);
  if (fst->Start() == kNoStateId) return;

  // New super final state
  const StateId sFinal = fst->AddState();

  std::vector<Arc> new_arcs_from_init;
  for (StateIterator<Fst> sit(*fst); !sit.Done(); sit.Next()) {
    const StateId u = sit.Value();
    if (u == sFinal || u == fst->Start()) continue;
    KALDI_ASSERT(u < state_group.size());
    const auto gu = state_group[u];
    std::vector<Arc> new_arcs_from_u;  // New arcs from state u

    // Add arc to the new (unique) final state, and make u not final
    if (fst->Final(u) != Weight::Zero()) {
      new_arcs_from_u.push_back(Arc(0, 0, fst->Final(u), sFinal));
      fst->SetFinal(u, Weight::Zero());
    }

    // Traverse current arcs from state u. When reaching a state v part of a
    // different group than u:
    //   - Remove current arc from the fst.
    //   - Add arc from the current state (u) to the final state, with cost =
    //     arc.weight * backward[v]. Since the state u is the final state of a
    //     group.
    //   - Add arc from the initial state to v, with cost =
    //     arc.weight * forward[u]. Since the state v is the start of a group.
    for (ArcIterator<Fst> ait(*fst, u); !ait.Done(); ait.Next()) {
      const Arc& arc = ait.Value();
      const StateId v = arc.nextstate;
      KALDI_ASSERT(v < state_group.size());
      const auto gv = state_group[v];
      if (gu != gv) {
        new_arcs_from_u.push_back(
            Arc(0, 0, fst::Times(arc.weight, bw[v]), sFinal));
        fst->AddArc(
            fst->Start(),
            Arc(arc.ilabel, arc.olabel, fst::Times(arc.weight, fw[u]), v));
      } else {
        new_arcs_from_u.push_back(arc);
      }
    }

    // Delete all original arcs from u, and all the new arcs
    fst->DeleteArcs(u);
    for (const Arc& arc : new_arcs_from_u) {
      fst->AddArc(u, arc);
    }
  }

  // Final cost = -total_cost, so that paths are normalized in -log [0, 1]
  fst->SetFinal(sFinal, Weight::One());
  // Remove epsilon symbols O(V^2 + V * E)
  RmEpsilon(fst);
  // Remove unnecessary states/arcs
  Connect(fst);
  // Push weights to the toward the initial state. This speeds up n-best list
  // retrieval.
  Push<Arc>(fst, fst::REWEIGHT_TO_INITIAL, fst::kPushWeights);
}

}  // namespace fst

#endif  // KALDI_FSTEXT_FSTEXT_UTILS2_H_

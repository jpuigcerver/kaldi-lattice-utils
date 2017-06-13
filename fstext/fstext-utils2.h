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
    const Fst<Arc>& ifst, MutableFst<Arc>* ofst, bool use_input = false,
    std::vector<size_t>* state_input_length = nullptr) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  if (ifst.Properties(kAcyclic, true) != kAcyclic) {
    KALDI_ERR << "DisambiguateStateInputSequenceLength() only supports acyclic "
              << "transducers";
  }

  ofst->DeleteStates();
  if (state_input_length != nullptr) state_input_length->clear();
  if (ifst.Start() < 0) return 0;


  // In the new fst, the states are represented by tuples:
  // (state, input length) from the old fst.
  // This map is used to map these tuples to StateId in the output fst.
  std::map<std::tuple<StateId, size_t>, StateId> state_map;
  // This queue is used to build the output fst.
  std::queue<std::tuple<StateId, StateId, size_t>> Q;

  // Initial state.
  ofst->SetStart(ofst->AddState());
  state_map[std::make_tuple(ifst.Start(), 0)] = 0;
  Q.push(std::make_tuple(ifst.Start(), ofst->Start(), 0));
  size_t max_len = 0;

  while (!Q.empty()) {
    const StateId u  = std::get<0>(Q.front());  // StateId in the old fst
    const StateId u2 = std::get<1>(Q.front());  // StateId in the new fst
    const size_t len = std::get<2>(Q.front());  // Input length to new state
    Q.pop();

    // Update the maximum length seen so far.
    if (max_len < len) max_len = len;

    // NOTE: Any state in the new fst is only added once the queue, and they
    // are added in the FIFO in the same order as they are created, Thus, we
    // can just push_back() the sequence length to the state_input_length
    // vector.
    if (state_input_length != nullptr) { state_input_length->push_back(len); }

    // Traverse edges from state u
    for (ArcIterator< Fst<Arc> > aiter(ifst, u); !aiter.Done(); aiter.Next()) {
      auto arc = aiter.Value();
      // If the arc is labeled with an epsilon, then the length of the sequence
      // is no incresead. Otherwise, it is increased.
      const size_t next_len =
          ((use_input ? arc.ilabel : arc.olabel) == 0) ? (len) : (len + 1);
      const auto state_tuple = std::make_tuple(arc.nextstate, next_len);
      const auto state_it = state_map.find(state_tuple);
      if (state_it == state_map.end()) {
        // If the tuple was never seen before, create a new state in the output
        // fst and add it to the map.
        const StateId nextstate2 = ofst->AddState();
        Q.push(std::make_tuple(arc.nextstate, nextstate2, next_len));
        state_map[state_tuple] = arc.nextstate = nextstate2;
      } else {
        // Otherwise, just use the previously assigned StateId.
        arc.nextstate = state_it->second;
      }
      // Add arc to the output fst
      ofst->AddArc(u2, arc);
    }

    // Set final weight for the new state in the output fst.
    ofst->SetFinal(u2, ifst.Final(u));
  }

  return max_len;
}

}  // namespace fst

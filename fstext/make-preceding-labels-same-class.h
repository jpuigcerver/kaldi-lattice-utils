#ifndef KALDI_LATTICE_UTILS_FSTEXT_MAKE_PRECEDING_LABELS_SAME_CLASS_H_
#define KALDI_LATTICE_UTILS_FSTEXT_MAKE_PRECEDING_LABELS_SAME_CLASS_H_

#include <queue>
#include <tuple>
#include <unordered_map>

#include "base/kaldi-error.h"
#include "fst/fstlib.h"
#include "util/tuple-hash.h"

namespace fst {

struct PrecedingLabelsSameClassOptions {
  bool use_input;
  bool propagate_epsilon_class;
  PrecedingLabelsSameClassOptions()
      : use_input(false), propagate_epsilon_class(false) {}
};

template <typename Arc, typename ClassType, typename F>
void MakePrecedingLabelsSameClass(
    const F& f,
    const Fst<Arc>& ifst,
    MutableFst<Arc>* ofst,
    std::vector<ClassType>* state_class,
    const PrecedingLabelsSameClassOptions& opts =
    PrecedingLabelsSameClassOptions()) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;

  KALDI_ASSERT(ofst != nullptr);
  KALDI_ASSERT(state_class != nullptr);

  ofst->DeleteStates();
  state_class->clear();
  if (ifst.Start() == kNoStateId) {
    return;
  }

  std::unordered_map<
      std::tuple<StateId, ClassType>,
      StateId,
      kaldi::hash<std::tuple<StateId, ClassType>>
  > state_map;
  std::queue< std::tuple<StateId, ClassType> > Q;
  const auto c_eps = f(0);
  Q.emplace(ifst.Start(), c_eps);

  // Create the initial state of the output FST, with class c_eps
  // (i.e. the same class as the epsilon label)
  ofst->SetStart(ofst->AddState());
  state_map.emplace(make_tuple(ifst.Start(), c_eps), ofst->Start());
  state_class->resize(1, c_eps);

  while (!Q.empty()) {
    // Current (StateId, ClassType) in the input FST
    const auto s = std::get<0>(Q.front());
    const auto c = std::get<1>(Q.front());
    // StateId in the output FST (this always exists)
    const auto v = state_map.find(Q.front())->second;
    Q.pop();

    // Set the final weight of the output state
    if (ifst.Final(s) != Weight::Zero()) {
      ofst->SetFinal(v, ifst.Final(s));
    }

    // Traverse arcs from the input FST
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      Arc arc = ait.Value();
      // Class of the arc's label
      const auto c_arc = f(opts.use_input ? arc.ilabel : arc.olabel);
      // Pair (StateId, Class) of the destination state
      const auto c_nstate =
          opts.propagate_epsilon_class ? (c_arc != c_eps ? c_arc : c) : c_arc;
      const auto t = std::make_tuple(arc.nextstate, c_nstate);
      auto it = state_map.find(t);
      if (it == state_map.end()) {
        // Destination (StateId, Class) is not represented in the output FST.
        // Add a node in the output fst, set it's class and emplace the pair
        // to create new arcs from there.
        it = state_map.emplace_hint(it, t, ofst->AddState());
        state_class->push_back(c_nstate);
        Q.emplace(t);
      }

      // Add arc to the output FST.
      arc.nextstate = it->second;
      ofst->AddArc(v, arc);
    }
  }

  auto outprops = ifst.Properties(kFstProperties, false);
  ofst->SetProperties(outprops, kFstProperties);
}

template <typename Arc, typename ClassType, typename F>
void MakePrecedingLabelsSameClass(
    const F& f,
    MutableFst<Arc>* mfst,
    std::vector<ClassType>* state_class,
    const PrecedingLabelsSameClassOptions& opts =
    PrecedingLabelsSameClassOptions()) {
  VectorFst<Arc> tmp(*mfst);
  MakePrecedingLabelsSameClass(f, tmp, mfst, state_class, opts);
}

}  // namespace fst

#endif  // KALDI_LATTICE_UTILS_FSTEXT_MAKE_PRECEDING_LABELS_SAME_CLASS_H_

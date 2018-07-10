#ifndef KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_SAME_LABEL_CLASS_H_
#define KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_SAME_LABEL_CLASS_H_

#include <sstream>
#include <stack>
#include <unordered_set>
#include <vector>

#include "fst/fstlib.h"
#include "fst/symbol-table.h"
#include "fstext/fstext-utils2.h"
#include "fstext/make-preceding-symbols-same-class.h"
#include "util/tuple-hash.h"

namespace fst {

namespace internal {

// Get the key of a vector of labels from a symbol table, or add it and return
// the associated key.
template<typename Label>
Label GetOrAddLabelVectorToSymbolTable(const std::vector<Label> &labels,
                                       SymbolTable *table) {
  KALDI_ASSERT(table != nullptr);
  if (labels.empty()) return 0;

  const auto str = LabelVectorToString(labels, nullptr, kStringSeparator);
  if (table->Member(str)) {
    return table->Find(str);
  } else {
    return table->AddSymbol(str);
  }
}

}  // namespace internal

struct ExpandSubpathsOptions {
  size_t max_subpath_length;
  bool use_input;

  ExpandSubpathsOptions()
      : max_subpath_length(std::numeric_limits<size_t>::max()),
        use_input(false) {}

  ExpandSubpathsOptions(size_t max_length, bool use_input)
      : max_subpath_length(max_length), use_input(use_input) {}
};

template<typename Arc, typename ClassType>
size_t ExpandSubpathsWithSameLabelClass(
    const std::vector<ClassType> &state_class,
    const Fst<Arc> &ifst,
    MutableFst<Arc> *ofst,
    const std::set<ClassType> &non_expandable_classes = std::set<ClassType>(),
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef std::tuple<
  StateId,               // Initial state of the subpath
  StateId,               // Final state of the subpath
  Weight,                // Weight of the subpath
  std::vector<Label>,    // Input labels of the subpath
  std::vector<Label>     // Output labels of the subpath
  > SState;
  // Check requirements
  KALDI_ASSERT(ofst != nullptr);
  if (ifst.Properties(kAcyclic, true) != kAcyclic) {
    KALDI_ERR << "Input FST must be acyclic!";
  }

  ofst->DeleteStates();
  if (ifst.Start() == kNoStateId) return 0;

  // Output FST has at most as many states as the input FST
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = ofst->AddState();
  }
  ofst->SetStart(ifst.Start());
  const StateId super_final = ofst->AddState();
  ofst->SetFinal(super_final, Weight::One());

  // Get (or create) symbols tables for the output fsts
  SymbolTable *isyms = ofst->MutableInputSymbols();
  SymbolTable *osyms = ofst->MutableOutputSymbols();
  if (isyms == nullptr) { isyms = new SymbolTable; }
  if (osyms == nullptr) { osyms = new SymbolTable; }
  if (!isyms->Member(0)) { isyms->AddSymbol("0", 0); }
  if (!osyms->Member(0)) { osyms->AddSymbol("0", 0); }

  KALDI_ASSERT(state_class.size() > ifst.Start());
  const auto c_eps = state_class[ifst.Start()];

  std::stack<SState> S;
  S.emplace(ifst.Start(),
            ifst.Start(),
            Weight::One(),
            std::vector < Label > {},
            std::vector < Label > {});

  std::unordered_set <
      std::tuple < StateId, size_t >,
      kaldi::hash < std::tuple < StateId, size_t >>
      > expanded_from_arc;

  // Convert subpaths to individual arcs
  size_t num_arcs = 0;
  while (!S.empty()) {
    const auto s0 = std::get<0>(S.top());  // initial state in the subpath
    const auto s1 = std::get<1>(S.top());  // current state in the subpath
    const auto w = std::get<2>(S.top());   // weight of the subpath
    auto ilbls = std::get<3>(S.top());     // input labels in the subpath
    auto olbls = std::get<4>(S.top());     // output labels in the subpath
    S.pop();

    KALDI_ASSERT(s1 < state_class.size());
    const auto cp = state_class[s1]; // class of the subpath

    bool add_arc = false;
    for (ArcIterator<Fst<Arc>> ait(ifst, s1); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto &label = opts.use_input ? arc.ilabel : arc.olabel;
      KALDI_ASSERT(arc.nextstate < state_class.size());
      const auto &cns =
          state_class[arc.nextstate] != c_eps
          ? state_class[arc.nextstate] : cp;
      if (cp == c_eps ||
          (cp == cns && non_expandable_classes.count(cns) == 0)) {
        if (arc.ilabel) ilbls.push_back(arc.ilabel);
        if (arc.olabel) olbls.push_back(arc.olabel);
        const auto curr_length = opts.use_input ? ilbls.size() : olbls.size();
        if (curr_length < opts.max_subpath_length) {
          S.emplace(s0, arc.nextstate, Times(w, arc.weight), ilbls, olbls);
        }
        if (arc.ilabel) ilbls.pop_back();
        if (arc.olabel) olbls.pop_back();
      } else {
        add_arc = true;
        if (!expanded_from_arc.count(make_tuple(s1, ait.Position()))) {
          std::vector<Label> new_ilbls, new_olbls;
          if (arc.ilabel) new_ilbls.push_back(arc.ilabel);
          if (arc.olabel) new_olbls.push_back(arc.olabel);
          const auto curr_length =
              opts.use_input ? new_ilbls.size() : new_olbls.size();
          if (curr_length < opts.max_subpath_length) {
            S.emplace(s1, arc.nextstate, arc.weight, new_ilbls, new_olbls);
            expanded_from_arc.emplace(s1, ait.Position());
          }
        }
      }
    }

    if (s0 != s1 && add_arc) {
      auto ilabel = internal::GetOrAddLabelVectorToSymbolTable(ilbls, isyms);
      auto olabel = internal::GetOrAddLabelVectorToSymbolTable(olbls, osyms);
      ofst->AddArc(s0, Arc(ilabel, olabel, w, s1));
      ++num_arcs;
    }

    if (ifst.Final(s1) != Weight::Zero()) {
      auto ilabel =
          internal::GetOrAddLabelVectorToSymbolTable(ilbls, isyms);
      auto olabel =
          internal::GetOrAddLabelVectorToSymbolTable(olbls, osyms);
      ofst->AddArc(
          s0, Arc(ilabel, olabel, Times(w, ifst.Final(s1)), super_final));
      ++num_arcs;
    }
  }

  ofst->SetInputSymbols(isyms);
  ofst->SetOutputSymbols(osyms);
  Connect(ofst);
  return num_arcs;
}

template<typename Arc, typename ClassType, typename F>
size_t ExpandSubpathsWithSameLabelClass(
    const F &f,
    const Fst<Arc> &ifst,
    MutableFst<Arc> *ofst,
    const std::set<ClassType> &non_expandable_classes = std::set<ClassType>(),
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  VectorFst<Arc> tmp;
  std::vector<ClassType> state_class;
  PrecedingSymbolsSameClassOptions opts2;
  opts2.use_input = opts.use_input;
  opts2.propagate_epsilon_class = true;
  MakePrecedingSymbolsSameClass<Arc, ClassType, F>(
      f, ifst, &tmp, &state_class, opts2);
  return ExpandSubpathsWithSameLabelClass(
      state_class, tmp, ofst, non_expandable_classes, opts);
}

template<typename Arc, typename ClassType, typename F>
size_t ExpandSubpathsWithSameLabelClass(
    const F &f,
    MutableFst<Arc> *mfst,
    const std::set<ClassType> &non_expandable_classes = std::set<ClassType>(),
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  VectorFst<Arc> tmp;
  std::vector<ClassType> state_class;
  PrecedingSymbolsSameClassOptions opts2;
  opts2.use_input = opts.use_input;
  opts2.propagate_epsilon_class = true;
  MakePrecedingSymbolsSameClass<Arc, ClassType, F>(
      f, *mfst, &tmp, &state_class, opts2);
  return ExpandSubpathsWithSameLabelClass<Arc, ClassType>(
      state_class, tmp, mfst, non_expandable_classes, opts);
}

}  // namespace fst

#endif  // KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_SAME_LABEL_CLASS_H_

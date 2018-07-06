#ifndef KALDI_LATTICE_UTILS_EXPAND_SUBPATHS_FST_H_
#define KALDI_LATTICE_UTILS_EXPAND_SUBPATHS_FST_H_

#include <sstream>
#include <stack>
#include <vector>

#include <fst/fstlib.h>

#include "fstext/fstext-utils2.h"
#include "fstext/label_group.h"

namespace fst {

namespace internal {

// Compute |S1 - S2|
template<typename Label>
size_t CountSetDifference(const std::set<Label> &set1,
                          const std::set<Label> &set2) {
  size_t diff = set1.size();
  for (Label label : set2) {
    if (set1.count(label)) {
      --diff;
      if (diff == 0) {
        break;
      }
    }
  }
  return diff;
}

// Convert a vector of labels to a string.
// e.g.: 1 2 3 9 12 -> "1_2_3_9_12"
template<typename Label>
std::string LabelVectorToString(const std::vector<Label> &labels) {
  std::string str;
  if (!labels.empty()) {
    str = std::to_string(labels[0]);
    for (size_t i = 1; i < labels.size(); ++i) {
      str += kStringSeparator + std::to_string(labels[i]);
    }
  }
  return str;
}

// Get the key of a vector of labels from a symbol table, or add it and return
// the associated key.
template<typename Label>
Label GetOrAddLabelVectorToSymbolTable(const std::vector<Label> &labels,
                                       SymbolTable *stable) {
  KALDI_ASSERT(stable != nullptr);
  if (labels.empty()) return 0;

  const auto str = LabelVectorToString(labels);
  auto key = stable->Find(str);
  if (key == kNoSymbol) {
    return stable->AddSymbol(str);
  } else {
    return key;
  }
}

}  // namespace internal

struct ExpandSubpathsOptions {
  size_t max_subpath_length;
  bool use_inputs;

  ExpandSubpathsOptions()
      : max_subpath_length(std::numeric_limits<size_t>::max()),
        use_inputs(false) {}

  ExpandSubpathsOptions(size_t max_length, bool use_inputs)
      : max_subpath_length(max_length), use_inputs(use_inputs) {}
};

template<typename Arc>
size_t ExpandSubpathsByLabelGroup(
    const Fst<Arc> &ifst, const LabelGroup<typename Arc::Label> &label_group,
    MutableFst<Arc> *ofst, const ExpandSubpathsOptions &opts) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef std::tuple<
  StateId,               // Initial state of the subpath
  StateId,               // Final state of the subpath
  Label,                 // Label of the states in the subpath
  Weight,                // Weight of the subpath
  std::vector<Label>,    // Input labels of the subpath
  std::vector<Label>     // Output labels of the subpath
  > SState;
  KALDI_ASSERT(ofst != nullptr);
  if (!ifst.Properties(kTopSorted, true)) {
    KALDI_ERR << "Input FST must be topologically sorted.";
  }

  // Output FST has at most as many states as the input FST
  ofst->DeleteStates();
  std::vector<std::set<Label>> state_input_groups{{0}};
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    // Add corresponding state to the output fst, including final weight
    const auto s = ofst->AddState();
    if (ifst.Final(s) != Weight::Zero()) {
      ofst->SetFinal(s, ifst.Final(s));
    }

    // Keep track of the different label groups entering each state.
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto &label = opts.use_inputs ? arc.ilabel : arc.olabel;
      if (arc.nextstate >= state_input_groups.size()) {
        state_input_groups.resize(arc.nextstate + 1);
      }
      if (label == 0) {
        // Propagate input groups from the previous state
        state_input_groups[arc.nextstate].insert(state_input_groups[s].begin(),
                                                 state_input_groups[s].end());
      } else {
        state_input_groups[arc.nextstate].insert(label_group[label]);
      }
    }
  }
  if (ofst->NumStates() == 0) {
    return 0;
  }
  ofst->SetStart(ifst.Start());

  // Get (or create) symbols tables for the output fsts
  SymbolTable *isyms = ofst->MutableInputSymbols();
  SymbolTable *osyms = ofst->MutableOutputSymbols();
  if (isyms == nullptr) { isyms = new SymbolTable; }
  if (osyms == nullptr) { osyms = new SymbolTable; }
  if (!isyms->Member(0)) { isyms->AddSymbol("<eps>", 0); }
  if (!osyms->Member(0)) { osyms->AddSymbol("<eps>", 0); }

  std::stack<SState> S;
  // Epsilon group always starts from the start state.
  S.emplace(ifst.Start(),
            ifst.Start(),
            0,
            Weight::One(),
            std::vector < Label > {},
            std::vector < Label > {});

  // Arcs in the output fst start from any arc which connects two states of
  // different groups.
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto label = opts.use_inputs ? arc.ilabel : arc.olabel;
      const auto group = label_group[label];
      if (group != 0
          && internal::CountSetDifference(state_input_groups[s], {0, group})) {
        std::vector<Label> ilabels, olabels;
        if (arc.ilabel != 0) ilabels.push_back(arc.ilabel);
        if (arc.olabel != 0) olabels.push_back(arc.olabel);
        S.emplace(s, arc.nextstate, group, arc.weight, ilabels, olabels);
      }
    }
  }

  // Convert subpaths to individual arcs
  size_t num_arcs = 0;
  while (!S.empty()) {
    const auto s0 = std::get<0>(S.top());  // initial state in the subpath
    const auto s1 = std::get<1>(S.top());  // current state in the subpath
    const auto path_g = std::get<2>(S.top());  // group of the subpath
    const auto w = std::get<3>(S.top());   // weight of the subpath
    auto ilbls = std::get<4>(S.top());     // input labels in the subpath
    auto olbls = std::get<5>(S.top());     // output labels in the subpath
    S.pop();

    bool add_arc = false;
    for (ArcIterator<Fst<Arc>> ait(ifst, s1); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto &label = opts.use_inputs ? arc.ilabel : arc.olabel;
      const auto &arc_g = label_group[label];
      if (path_g == 0 || arc_g == 0 || path_g == arc_g) {
        if (arc.ilabel) ilbls.push_back(arc.ilabel);
        if (arc.olabel) olbls.push_back(arc.olabel);
        const auto curr_length = opts.use_inputs ? ilbls.size() : olbls.size();
        if (curr_length < opts.max_subpath_length) {
          S.emplace(s0,
                    arc.nextstate,
                    path_g == 0 ? arc_g : path_g,
                    Times(w, arc.weight),
                    ilbls,
                    olbls);
        }
        if (arc.ilabel) ilbls.pop_back();
        if (arc.olabel) olbls.pop_back();
      } else {
        add_arc = true;
      }
    }

    if (s0 != s1 && (add_arc || ofst->Final(s1) != Weight::Zero())) {
      auto ilabel = internal::GetOrAddLabelVectorToSymbolTable(ilbls, isyms);
      auto olabel = internal::GetOrAddLabelVectorToSymbolTable(olbls, osyms);
      ofst->AddArc(s0, Arc(ilabel, olabel, w, s1));
      ++num_arcs;
    }
  }

  ofst->SetInputSymbols(isyms);
  ofst->SetOutputSymbols(osyms);
  Connect(ofst);
  return num_arcs;
}

template<typename Arc>
size_t ExpandSubpathsByStateGroup(
    const Fst<Arc> &ifst, const std::vector<typename Arc::Label> &state_group,
    MutableFst<Arc> *ofst, const ExpandSubpathsOptions &opts) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef std::tuple<
  StateId,               // Initial state of the subpath
  StateId,               // Final state of the subpath
  Label,                 // Label of the states in the subpath
  Weight,                // Weight of the subpath
  std::vector<Label>,    // Input labels of the subpath
  std::vector<Label>     // Output labels of the subpath
  > SState;
  KALDI_ASSERT(ofst != nullptr);
  if (!ifst.Properties(kTopSorted, true)) {
    KALDI_ERR << "Input FST must be topologically sorted.";
  }

  // Output FST has at most as many states as the input FST
  ofst->DeleteStates();
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    // Add corresponding state to the output fst, including final weight
    const auto s = ofst->AddState();
    if (ifst.Final(s) != Weight::Zero()) {
      ofst->SetFinal(s, ifst.Final(s));
    }
  }
  if (ofst->NumStates() == 0) {
    return 0;
  }
  ofst->SetStart(ifst.Start());

  // Get (or create) symbols tables for the output fsts
  SymbolTable *isyms = ofst->MutableInputSymbols();
  SymbolTable *osyms = ofst->MutableOutputSymbols();
  if (isyms == nullptr) { isyms = new SymbolTable; }
  if (osyms == nullptr) { osyms = new SymbolTable; }
  if (!isyms->Member(0)) { isyms->AddSymbol("<eps>", 0); }
  if (!osyms->Member(0)) { osyms->AddSymbol("<eps>", 0); }

  std::stack<SState> S;
  // Epsilon group always starts from the start state.
  S.emplace(ifst.Start(),
            ifst.Start(),
            0,
            Weight::One(),
            std::vector < Label > {},
            std::vector < Label > {});


  //////////////////////////////////////////////////////////
  //////////////////////// TO DO ///////////////////////////
  //////////////////////////////////////////////////////////

  // Arcs in the output fst start from any arc which connects two states of
  // different groups.
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto label = opts.use_inputs ? arc.ilabel : arc.olabel;
      const auto group = label_group[label];
      if (group != 0
          && internal::CountSetDifference(state_input_groups[s], {0, group})) {
        std::vector<Label> ilabels, olabels;
        if (arc.ilabel != 0) ilabels.push_back(arc.ilabel);
        if (arc.olabel != 0) olabels.push_back(arc.olabel);
        S.emplace(s, arc.nextstate, group, arc.weight, ilabels, olabels);
      }
    }
  }

  // Convert subpaths to individual arcs
  size_t num_arcs = 0;
  while (!S.empty()) {
    const auto s0 = std::get<0>(S.top());  // initial state in the subpath
    const auto s1 = std::get<1>(S.top());  // current state in the subpath
    const auto path_g = std::get<2>(S.top());  // group of the subpath
    const auto w = std::get<3>(S.top());   // weight of the subpath
    auto ilbls = std::get<4>(S.top());     // input labels in the subpath
    auto olbls = std::get<5>(S.top());     // output labels in the subpath
    S.pop();

    bool add_arc = false;
    for (ArcIterator<Fst<Arc>> ait(ifst, s1); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto &label = opts.use_inputs ? arc.ilabel : arc.olabel;
      const auto &arc_g = label_group[label];
      if (path_g == 0 || arc_g == 0 || path_g == arc_g) {
        if (arc.ilabel) ilbls.push_back(arc.ilabel);
        if (arc.olabel) olbls.push_back(arc.olabel);
        const auto curr_length = opts.use_inputs ? ilbls.size() : olbls.size();
        if (curr_length < opts.max_subpath_length) {
          S.emplace(s0,
                    arc.nextstate,
                    path_g == 0 ? arc_g : path_g,
                    Times(w, arc.weight),
                    ilbls,
                    olbls);
        }
        if (arc.ilabel) ilbls.pop_back();
        if (arc.olabel) olbls.pop_back();
      } else {
        add_arc = true;
      }
    }

    if (s0 != s1 && (add_arc || ofst->Final(s1) != Weight::Zero())) {
      auto ilabel = internal::GetOrAddLabelVectorToSymbolTable(ilbls, isyms);
      auto olabel = internal::GetOrAddLabelVectorToSymbolTable(olbls, osyms);
      ofst->AddArc(s0, Arc(ilabel, olabel, w, s1));
      ++num_arcs;
    }
  }

  ofst->SetInputSymbols(isyms);
  ofst->SetOutputSymbols(osyms);
  Connect(ofst);
  return num_arcs;
}

}  // namespace fst

#endif  // KALDI_LATTICE_UTILS_EXPAND_SUBPATHS_FST_H_

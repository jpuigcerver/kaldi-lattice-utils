#ifndef KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_BETWEEN_DELIMITERS_H_
#define KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_BETWEEN_DELIMITERS_H_

#include "fstext/expand-subpaths-labels-same-class.h"

namespace fst {

namespace  {

struct StateInfo {
  bool has_inp_delim;
  bool has_inp_regular;
  bool has_out_delim;
  bool has_out_regular;

  StateInfo()
      : has_inp_delim(false), has_inp_regular(false),
        has_out_delim(false), has_out_regular(false) {}

  StateInfo(bool has_inp_delim, bool has_inp_regular,
            bool has_out_delim, bool has_out_regular)
      : has_inp_delim(has_inp_delim), has_inp_regular(has_inp_regular),
        has_out_delim(has_out_delim), has_out_regular(has_out_regular) {}
};

template <typename Arc>
bool CanUseExpandSubpathsBetweenDelimitersSpecial(
    const std::set<typename Arc::Label>& delimiters,
    const Fst<Arc>& ifst,
    const bool use_input) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  std::unordered_map<StateId, StateInfo> state_info;
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();

    bool has_out_delim = false, has_out_regular = false;
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto& arc = ait.Value();
      const auto label = use_input ? arc.ilabel : arc.olabel;
      if (label == 0) return false;
      const auto is_delim = (delimiters.count(label) != 0);

      // update input info for state arc.nextstate
      auto it = state_info.find(arc.nextstate);
      if (it == state_info.end()) {
        state_info.emplace_hint(
            state_info.end(),
            arc.nextstate,
            StateInfo{is_delim, !is_delim, false, false});
      } else {
        it->second.has_inp_delim = it->second.has_inp_delim || is_delim;
        it->second.has_inp_regular = it->second.has_inp_regular || !is_delim;
      }

      // update output info for state s
      if (is_delim) { has_out_delim = true; } else { has_out_regular = true; }
    }

    // write output info for state s
    auto it = state_info.find(s);
    if (it == state_info.end()) {
      state_info.emplace_hint(
          state_info.end(),
          s,
          StateInfo{false, false, has_out_delim, has_out_regular});
    } else {
      it->second.has_out_delim = has_out_delim;
      it->second.has_out_regular = has_out_regular;
    }
  }

  for (const auto& kv : state_info) {
    if (kv.second.has_inp_regular && kv.second.has_inp_delim &&
        kv.second.has_out_regular &&
        (kv.second.has_out_delim || ifst.Final(kv.first) != Weight::Zero())) {
      return false;
    }
  }

  return true;
}

template<typename Arc>
void ExpandSubpathsBetweenDelimitersSpecial(
    const std::set<typename Arc::Label> &delimiters,
    const Fst <Arc> &ifst,
    MutableFst <Arc> *ofst,
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(ofst != nullptr);
  if (ifst.Properties(kAcyclic, true) != kAcyclic) {
    KALDI_ERR << "Input FST must be acyclic!";
  }

  ofst->DeleteStates();
  if (ifst.Start() == kNoStateId) return;

  // Get (or create) symbols tables for the output fsts
  SymbolTable *isyms = ofst->MutableInputSymbols();
  SymbolTable *osyms = ofst->MutableOutputSymbols();
  if (isyms == nullptr) { isyms = new SymbolTable; }
  if (osyms == nullptr) { osyms = new SymbolTable; }
  if (!isyms->Member(0)) { isyms->AddSymbol("0", 0); }
  if (!osyms->Member(0)) { osyms->AddSymbol("0", 0); }

  // Output fst has, at most, as many states as the input fst.
  ofst->DeleteStates();
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    ofst->AddState();
  }
  ofst->SetStart(ifst.Start());

  // Add start state to the stack, new words start from these states.
  std::unordered_set<StateId> add_to_initial_stack;
  add_to_initial_stack.insert(ifst.Start());

  // Traverse all arcs from the fst to detect where new words may start.
  // We will add these states to the stack, and keep the arcs with any of the
  // word delimiters, in order to keep these as separate words.
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const StateId s = sit.Value();
    ofst->SetFinal(s, ifst.Final(s));

    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const Arc &arc = ait.Value();
      const Label &label = opts.use_input ? arc.ilabel : arc.olabel;
      // State arc.nextstate is the start of a new word
      if (delimiters.count(label) != 0) {
        // Get ilabel and olabel corresponding to this arc in the output fst
        std::vector<Label> ilabels, olabels;
        if (arc.ilabel != 0) ilabels.push_back(arc.ilabel);
        if (arc.olabel != 0) olabels.push_back(arc.olabel);
        const auto ilabel = GetOrAddLabelVectorToSymbolTable(ilabels, isyms);
        const auto olabel = GetOrAddLabelVectorToSymbolTable(olabels, osyms);
        // Add arc to the output fst and add it to the stack
        ofst->AddArc(s, Arc(ilabel, olabel, arc.weight, arc.nextstate));
        add_to_initial_stack.insert(arc.nextstate);
      }
    }
  }

  // Initialize stack
  std::stack< std::tuple<StateId, StateId, Path<Label, Weight>> > S;
  for (StateId s : add_to_initial_stack) {
    S.emplace(s, s, Path<Label, Weight>());
  }

  while (!S.empty()) {
    const auto i = std::get<0>(S.top());
    const auto j = std::get<1>(S.top());
    const auto p = std::get<2>(S.top());
    S.pop();

    bool add_arc = false;
    for (ArcIterator<Fst<Arc>> ait(ifst, j); !ait.Done(); ait.Next()) {
      const Arc& arc = ait.Value();
      const Label label = opts.use_input ? arc.ilabel : arc.olabel;
      if (delimiters.count(label) == 0) {
        // Continue expanding path from (i, arc.nextstate)
        const auto new_p =
            Path<Label, Weight>(p, arc.weight, arc.ilabel, arc.olabel);
        if (new_p.Length(opts.use_input) <= opts.max_subpath_length) {
          S.emplace(i, arc.nextstate, new_p);
        }
      } else {
        add_arc = true;
      }
    }

    // If j is final or has any output arc with a delimiter symbol,
    // add arc from i to j to the output fst, with the cumulative weight
    // and input/output labels.
    if (i != j && (add_arc || ifst.Final(j) != Weight::Zero())) {
      const auto ilabel = GetOrAddLabelVectorToSymbolTable(p.ilabels, isyms);
      const auto olabel = GetOrAddLabelVectorToSymbolTable(p.olabels, osyms);
      ofst->AddArc(i, Arc(ilabel, olabel, p.weight, j));
    }
  }

  // Remove unnecessary states and arcs
  Connect(ofst);
  ofst->SetInputSymbols(isyms);
  ofst->SetOutputSymbols(osyms);
}

}  // namespace

template<typename Arc>
void ExpandSubpathsBetweenDelimiters(
    const std::set<typename Arc::Label> &delimiters,
    const Fst<Arc>& ifst,
    MutableFst <Arc> *ofst,
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  typedef typename Arc::Label Label;

  if (CanUseExpandSubpathsBetweenDelimitersSpecial(
      delimiters, ifst, opts.use_input)) {
    KALDI_VLOG(4) << "Using special subpath expansion between delimiters...";
    ExpandSubpathsBetweenDelimitersSpecial(delimiters, ifst, ofst, opts);
  } else {
    KALDI_VLOG(4) << "Using general subpath expansion between delimiters...";
    auto f = [&delimiters](Label label) -> int {
      return label == 0 ? 0 : (delimiters.count(label) == 0 ? 1 : 2);
    };
    ExpandSubpathsLabelsSameClass<Arc, int>(
        f, ifst, ofst, std::set<int>{2}, opts);
  }
}

template<typename Arc>
void ExpandSubpathsBetweenDelimiters(
    const std::set<typename Arc::Label> &delimiters,
    MutableFst <Arc> *mfst,
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  VectorFst<Arc> tmp(*mfst);
  ExpandSubpathsBetweenDelimiters(delimiters, tmp, mfst, opts);
}

}  // namespace fst

#endif // KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_BETWEEN_DELIMITERS_H_

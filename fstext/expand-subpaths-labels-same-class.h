#ifndef KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_SAME_LABEL_CLASS_H_
#define KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_SAME_LABEL_CLASS_H_

#include <sstream>
#include <stack>
#include <unordered_set>
#include <vector>

#include "fst/fstlib.h"
#include "fst/symbol-table.h"
#include "fstext/fstext-utils2.h"
#include "util/tuple-hash.h"

namespace fst {

struct ExpandSubpathsOptions {
  size_t max_subpath_length;
  bool use_input;

  ExpandSubpathsOptions()
      : max_subpath_length(std::numeric_limits<size_t>::max()),
        use_input(false) {}

  ExpandSubpathsOptions(size_t max_length, bool use_input)
      : max_subpath_length(max_length), use_input(use_input) {}
};

namespace {

// Get the key of a vector of labels from a symbol table, or add it and return
// the associated key.
template<typename Label>
Label GetOrAddLabelVectorToSymbolTable(const std::vector <Label> &labels,
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

template<typename Label, typename Weight>
struct Path {
  Weight weight;
  std::vector <Label> ilabels;
  std::vector <Label> olabels;

  Path() : weight(Weight::One()), ilabels{}, olabels{} {}

  Path(const Path &other)
      : weight(other.weight), ilabels(other.ilabels), olabels(other.olabels) {}

  explicit Path(const Weight &w) : weight(w), ilabels{}, olabels{} {}

  Path(const Weight &w, const Label &ilabel, const Label &olabel)
      : weight(w), ilabels{}, olabels{} {
    if (ilabel != 0) ilabels.push_back(ilabel);
    if (olabel != 0) olabels.push_back(olabel);
  }

  Path(const Path<Label, Weight> &prev,
       const Weight &w,
       const Label &ilabel,
       const Label &olabel)
      : weight(Times(prev.weight, w)),
        ilabels(prev.ilabels),
        olabels(prev.olabels) {
    if (ilabel != 0) ilabels.push_back(ilabel);
    if (olabel != 0) olabels.push_back(olabel);
  }

  size_t Length(bool use_input) const {
    return use_input ? ilabels.size() : olabels.size();
  }
};

}  // namespace


template<typename Arc, typename ClassType, typename F>
void ExpandSubpathsLabelsSameClass(
    const F &f,
    const Fst <Arc> &ifst,
    MutableFst <Arc> *ofst,
    const std::set <ClassType> &non_expandable_classes = std::set<ClassType>(),
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
  ofst->SetStart(ofst->AddState());

  // Get (or create) symbols tables for the output fsts
  SymbolTable *isyms = ofst->MutableInputSymbols();
  SymbolTable *osyms = ofst->MutableOutputSymbols();
  if (isyms == nullptr) { isyms = new SymbolTable; }
  if (osyms == nullptr) { osyms = new SymbolTable; }
  if (!isyms->Member(0)) { isyms->AddSymbol("0", 0); }
  if (!osyms->Member(0)) { osyms->AddSymbol("0", 0); }

  const auto &c_eps = f(0);

  std::unordered_map <
      std::tuple < StateId, ClassType >,
      StateId,
      kaldi::hash < std::tuple < StateId, ClassType >> > M;
  M[std::make_tuple(ifst.Start(), c_eps)] = ofst->Start();
  std::vector <std::tuple<StateId, ClassType>> IM{
      std::make_tuple(ifst.Start(), c_eps)
  };

  std::stack <std::tuple<StateId, StateId, ClassType, Path<Label, Weight>>> S;
  S.emplace(ofst->Start(), ifst.Start(), c_eps, Path<Label, Weight>());

  std::unordered_set <
      std::tuple < StateId, size_t >,
      kaldi::hash < std::tuple < StateId, size_t >> > X;

  size_t num_arcs = 0;
  while (!S.empty()) {
    const auto i = std::get<0>(S.top());
    const auto j = std::get<1>(S.top());
    const auto c = std::get<2>(S.top());
    const auto p = std::get<3>(S.top());
    S.pop();

    bool add_arc = false;
    for (ArcIterator <Fst<Arc>> ait(ifst, j); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto label = opts.use_input ? arc.ilabel : arc.olabel;
      const auto c_arc = f(label) != c_eps ? f(label) : c;
      if (c == c_eps || (c == c_arc && !non_expandable_classes.count(c_arc))) {
        const auto new_p =
            Path<Label, Weight>(p, arc.weight, arc.ilabel, arc.olabel);
        if (new_p.Length(opts.use_input) <= opts.max_subpath_length) {
          S.emplace(i, arc.nextstate, c_arc, new_p);
        }
      } else {
        add_arc = true;
        const auto
            new_p = Path<Label, Weight>(arc.weight, arc.ilabel, arc.olabel);
        if (new_p.Length(opts.use_input) <= opts.max_subpath_length) {
          const auto t = std::make_tuple(j, c);
          auto it = M.find(t);
          if (it == M.end()) {
            it = M.emplace_hint(M.end(), t, ofst->AddState());
            IM.emplace_back(t);
          }
          if (X.find(std::make_tuple(it->second, ait.Position())) == X.end()) {
            S.emplace(it->second, arc.nextstate, c_arc, new_p);
            X.emplace(it->second, ait.Position());
          }
        }
      }
    }

    if (j != std::get<0>(IM[i])
        && (ifst.Final(j) != Weight::Zero() || add_arc)) {
      const auto t = std::make_tuple(j, c);
      auto it = M.find(t);
      if (it == M.end()) {
        it = M.emplace_hint(M.end(), t, ofst->AddState());
        IM.emplace_back(t);
      }
      auto ilabel =
          GetOrAddLabelVectorToSymbolTable(p.ilabels, isyms);
      auto olabel =
          GetOrAddLabelVectorToSymbolTable(p.olabels, osyms);
      ofst->AddArc(i, Arc(ilabel, olabel, p.weight, it->second));
      ++num_arcs;
    }
  }

  for (const auto &kv : M) {
    const auto &final_weight = ifst.Final(std::get<0>(kv.first));
    if (final_weight != Weight::Zero()) {
      ofst->SetFinal(kv.second, final_weight);
    }
  }

  Connect(ofst);
  ofst->SetInputSymbols(isyms);
  ofst->SetOutputSymbols(osyms);
}

template<typename Arc, typename ClassType, typename F>
void ExpandSubpathsLabelsSameClass(
    const F &f,
    MutableFst <Arc> *mfst,
    const std::set <ClassType> &non_expandable_classes = std::set<ClassType>(),
    const ExpandSubpathsOptions &opts = ExpandSubpathsOptions()) {
  VectorFst <Arc> tmp(*mfst);
  ExpandSubpathsLabelsSameClass<Arc, ClassType>(
      f, tmp, mfst, non_expandable_classes, opts);
}

}  // namespace fst

#endif  // KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_SAME_LABEL_CLASS_H_

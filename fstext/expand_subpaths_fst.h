#ifndef KALDI_LATTICE_UTILS_EXPAND_SUBPATHS_FST_H_
#define KALDI_LATTICE_UTILS_EXPAND_SUBPATHS_FST_H_

#include <stack>
#include <vector>

#include <fst/fstlib.h>

#include "fstext/fstext-utils2.h"

namespace fst {

namespace internal {

template<typename Label>
class VectorLabelHash {
  size_t operator()(const std::vector <Label> &v) const {
    hash <Label> hashe;
    size_t h = 0;
    for (const auto &s : v) {
      h ^= hashe(s) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

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
void ExpandSubpathsByStateGroup(
    const Fst <Arc> &ifst, const std::vector<typename Arc::Label> &state_group,
    MutableFst <Arc> *ofst,
    std::vector <std::vector<typename Arc::Label>> *ilabel_map,
    std::vector <std::vector<typename Arc::Label>> *olabel_map,
    const ExpandSubpathsOptions &opts) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef std::tuple <
  StateId,               // Initial state of the subpath
  StateId,               // Final state of the subpath
  Label,                 // Label of the states in the subpath
  Weight,                // Weight of the subpath
  std::vector<Label>,    // Input labels of the subpath
  std::vector<Label>     // Output labels of the subpath
  > SState;
  KALDI_ASSERT(ofst != nullptr);
  KALDI_ASSERT(ilabel_map != nullptr);
  KALDI_ASSERT(olabel_map != nullptr);

  ilabel_map->clear();
  olabel_map->clear();

  // Output FST has at most as many states as the input FST
  ofst->DeleteStates();
  for (StateIterator <Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = ofst->AddState();
    if (ifst.Final(s) != Weight::Zero()) {
      ofst->SetFinal(s, ifst.Final(s));
    }
  }
  if (ofst->NumStates() == 0) {
    return;
  }
  ofst->SetStart(ifst.Start());

  std::stack <SState> S;
  // Epsilon group always starts from the start state.
  S.emplace(ifst.Start(), ifst.Start(), 0, Weight::One(), {}, {});

  // Arcs in the output fst start from any arc which connects two states of
  // different groups.
  for (StateIterator <Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    const auto gs = state_group[s];
    for (ArcIterator <Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      const auto gns = state_group[arc.nextstate];
      if (gs != 0 && gns != 0 && gs != gns) {
        std::vector <Label> ilabels, olabels;
        if (arc.ilabel != 0) ilabels.push_back(arc.ilabel);
        if (arc.olabel != 0) olabels.push_back(arc.olabel);
        S.emplace(s, arc.nextstate, gns, arc.weight, ilabels, olabels);
      }
    }
  }

  std::unorderd_map <std::vector<Label>, Label, VectorLabelHash> iseq2lbl,
      oseq2lbl;
  iseq2lbl[std::vector < Label > {}] = 0;
  oseq2lbl[std::vector < Label > {}] = 0;

  // Process subpaths
  while (!S.empty()) {
    const auto s0 = std::get<1>(S.top());  // initial state in the subpath
    const auto s1 = std::get<2>(S.top());  // current state in the subpath
    const auto g = std::get<3>(S.top());   // group of the subpath
    const auto w = std::get<4>(S.top());   // weight of the subpath
    auto ilbls = std::get<5>(S.top());     // input labels in the subpath
    auto olbls = std::get<6>(S.top());     // output labels in the subpath
    S.pop();

    bool add_arc = false;
    for (ArcIterator <Fst<Arc>> ait(ifst, ns); !ait.Done(); ait.Next()) {
      const auto &arc = ait.Value();
      if (gr == 0 ||
          state_group[arc.nextstate] == 0 ||
          gr == state_group[arc.nextstate]) {
        if (arc.ilabel) ilbls.push_back(arc.ilabel);
        if (arc.olabel) olbls.push_back(arc.olabel);
        const auto curr_length = opts.use_inputs ? ilbls.size() : olbls.size();
        if (curr_length < opts.max_subpath_length) {
          S.emplace(cs,
                    arc.nextstate,
                    gr == 0 ? state_group[arc.nextstate] : gr,
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

    if (s0 != s1 && (add_arc || ofst->Final(ns) != Weight::Zero())) {
      const auto new_ilabel =
          iseq2lbl.emplace_back(ilbls, iseq2lbl.size()).first->second;
      const auto new_olabel =
          oseq2lbl.emplace_back(olbls, oseq2lbl.size()).first->second;
      ofst->AddArc(s0, Arc(new_ilabel, new_olabel, w, s1));
    }
  }

  // Set mappings from Label to Vector<Label> (to recover original sequences
  // from label integers in the output fst)
  ilabel_map->resize(iseq2lbl.size());
  for (const auto& kv : iseq2lbl) {
    (*ilabel_map)[kv.second] = kv.first;
  }
  olabel_map->resize(oseq2lbl.size());
  for (const auto& kv : oseq2lbl) {
    (*olabel_map)[kv.second] = kv.first;
  }

  // Trim output FST to remove unreachable states
  Connect(ofst);
}

template<typename Arc>
void ExpandSubpathsByLabelGroup(
    const Fst <Arc> &ifst, const LabelGroup<typename Arc::Label>& label_group,
    MutableFst <Arc> *ofst,
    std::vector <std::vector<typename Arc::Label>> *ilabel_map,
    std::vector <std::vector<typename Arc::Label>> *olabel_map,
    const ExpandSubpathsOptions &opts) {
  std::vector<typename Arc::Label> state_group;
  DisambiguateStatesByInputLabelGroup(ifst, label_group, ofst, opts.use_input);
  ExpandSubpathsByStateGroup(ifst, state_group, ofst, ilabel_map, olabel_map, opts);
}

}  // namespace fst

#endif  // KALDI_LATTICE_UTILS_EXPAND_SUBPATHS_FST_H_

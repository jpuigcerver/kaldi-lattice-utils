#ifndef KALDI_LATTICE_UTILS_FSTEXT_FST_TEST_UTILS2_H_
#define KALDI_LATTICE_UTILS_FSTEXT_FST_TEST_UTILS2_H_

#include <algorithm>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "fst/fstlib.h"
#include "fstext/fstext-utils2.h"

namespace fst {

void GetBestPathsFromFst(
    const Fst<StdArc>& fst,
    std::vector<std::tuple<float, std::string, std::string>>* nbests,
    const int num_paths = std::numeric_limits<int>::max()) {
  VectorFst<StdArc> nbest_fst;
  ShortestPath(fst, &nbest_fst, num_paths);

  std::vector<VectorFst<StdArc>> vec_nbest_fst;
  ConvertNbestToVector(nbest_fst, &vec_nbest_fst);

  nbests->clear();
  for (const auto& nbest : vec_nbest_fst) {
    std::vector<StdArc::Label> isyms, osyms;
    StdArc::Weight cost;
    GetLinearSymbolSequence(nbest, &isyms, &osyms, &cost);
    nbests->emplace_back(
        cost.Value(),
        LabelVectorToString(isyms, fst.InputSymbols()),
        LabelVectorToString(osyms, fst.OutputSymbols()));
  }
  // Sort again, if two paths have the same cost, order w.r.t. input and
  // output labels
  std::sort(nbests->begin(), nbests->end());
}

bool AreBestPathsEqual(
    const std::vector<std::tuple<float, std::string, std::string>>& nbests1,
    const std::vector<std::tuple<float, std::string, std::string>>& nbests2,
    size_t* diff_pos = nullptr) {
  if(diff_pos) *diff_pos = 0;
  if (nbests1.size() != nbests2.size()) return false;

  for (size_t i = 0; i < std::min(nbests1.size(), nbests2.size()); ++i) {
    const auto& nb1 = nbests1[i];
    const auto& nb2 = nbests2[i];
    if (!kaldi::ApproxEqual(std::get<0>(nb1), std::get<0>(nb2)) ||
        std::get<1>(nb1) != std::get<1>(nb2) ||
        std::get<2>(nb1) != std::get<2>(nb2)) {
      if (diff_pos) *diff_pos = i + 1;
      return false;
    }
  }

  return true;
}

bool AreBestPathsEqual(
    const Fst<StdArc>& fst1, const Fst<StdArc>& fst2,
    size_t* diff_pos = nullptr,
    const int num_paths = std::numeric_limits<int>::max()) {
  std::vector<std::tuple<float, std::string, std::string>> nbests1, nbests2;
  GetBestPathsFromFst(fst1, &nbests1, num_paths);
  GetBestPathsFromFst(fst2, &nbests2, num_paths);
  return AreBestPathsEqual(nbests1, nbests2, diff_pos);
}


template <typename Arc>
void RenumberStatesFst(const Fst<Arc>& ifst, MutableFst<Arc>* ofst) {
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  ofst->DeleteStates();
  if (ifst.Start() == kNoStateId) return;

  std::unordered_map<StateId, StateId> state_map;
  ofst->SetStart(ofst->AddState());
  state_map[ifst.Start()] = ofst->Start();

  std::queue<StateId> Q;
  Q.push(ifst.Start());

  while (!Q.empty()) {
    const StateId u1 = Q.front();
    const StateId u2 = state_map.find(u1)->second;
    Q.pop();

    if (ifst.Final(u1) != Weight::Zero()) {
      ofst->SetFinal(u2, ifst.Final(u1));
    }

    for (ArcIterator<Fst<Arc>> ait(ifst, u1); !ait.Done(); ait.Next()) {
      Arc arc = ait.Value();
      auto it = state_map.find(arc.nextstate);
      if (it == state_map.end()) {
        it = state_map.emplace_hint(it, arc.nextstate, ofst->AddState());
        Q.push(arc.nextstate);
      }
      arc.nextstate = it->second;
      ofst->AddArc(u2, arc);
    }
  }
}

template <typename Arc>
void RenumberStatesFst(MutableFst<Arc>* mfst) {
  VectorFst<Arc> copy(*mfst);
  RenumberStatesFst(copy, mfst);
}

template <typename Label>
class RandomClass {
 public:
  RandomClass(const size_t num_classes = 8, const size_t num_labels = 1000)
      : label2class_(num_labels) {
    for (size_t i = 0; i < label2class_.size(); ++i) {
      label2class_[i] = (int)(rand() % num_classes) - 4;
    }
  }

  int operator()(Label label) const {
    return label2class_[label];
  }

  template <typename Arc>
  void PrintClasses(const Fst<Arc>& fst, bool use_input) const {
    std::set<Label> labels{0};
    for (StateIterator<Fst<Arc>> sit(fst); !sit.Done(); sit.Next()) {
      for (ArcIterator<Fst<Arc>> ait(fst, sit.Value()); !ait.Done(); ait.Next()) {
        labels.emplace(use_input ? ait.Value().ilabel : ait.Value().olabel);
      }
    }
    for (const auto& label : labels) {
      std::cerr << label << " -> " << (*this)(label) << std::endl;
    }
  }

 private:
  std::vector<int> label2class_;
};


}  // namespace fst

#endif  // KALDI_LATTICE_UTILS_FSTEXT_FST_TEST_UTILS2_H_

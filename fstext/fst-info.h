// fstext/fstext-info.h

// Copyright (c) 2018 Joan Puigcerver <joapuipe@upv.es>

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

#ifndef KALDI_LATTICE_UTILS_FST_INFO_H_
#define KALDI_LATTICE_UTILS_FST_INFO_H_

#include <iomanip>
#include <limits>

#include "fstext/fstext-lib.h"
#include "fst/script/info-impl.h"

namespace {
template<typename K>
void MapSetMax(std::map<K, int>& m, const K& k, int v) {
  auto it = m.find(k);
  if (it == m.end()) {
    m.emplace_hint(it, k, v);
  } else {
    it->second = std::max(it->second, v);
  }
}

template<typename K>
int MapGetDefault(const std::map<K, int>& m, const K& k, int v) {
  auto it = m.find(k);
  if (it == m.end()) return v;
  else return it->second;
}
}  // namespace

namespace fst {

template<typename FST>
long double ComputeNumberOfPaths(FST *fst) {
  typedef typename FST::Arc Arc;
  typedef typename Arc::Weight Weight;
  typedef typename FST::StateId StateId;

  if (fst->Properties(kAcyclic, true) & kAcyclic) {
    if (fst->Start() == kNoStateId) return 0;
    TopSort(fst);
    // (Log)number of paths arriving to each state.
    std::vector<long double> num_paths(fst->NumStates(), 0);
    num_paths[fst->Start()] = 1;
    // (Log)number of full paths.
    long double num_full_paths = 0;
    for (StateIterator<FST> siter(*fst); !siter.Done(); siter.Next()) {
      const StateId s = siter.Value();
      for (ArcIterator<FST> aiter(*fst, s); !aiter.Done();
           aiter.Next()) {
        const Arc &arc = aiter.Value();
        num_paths[arc.nextstate] += num_paths[s];
      }
      if (fst->Final(s) != Weight::Zero()) {
        num_full_paths += num_paths[s];
      }
    }
    return num_full_paths;
  } else {
    return std::numeric_limits<long double>::infinity();
  }
}

template <typename FST>
int ComputeMaxPathLength(FST* fst) {
  typedef typename FST::Arc Arc;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;

  if (fst->Properties(kAcyclic, true) & kAcyclic) {
    if (fst->Start() == kNoStateId) { return -1; }
    TopSort(fst);

    std::map<StateId, int> M;
    M[fst->Start()] = 0;
    for (StateIterator<FST> sit(*fst); !sit.Done(); sit.Next()) {
      const auto& s = sit.Value();
      const auto& l = M[s];
      for (ArcIterator<FST> ait(*fst, s); !ait.Done(); ait.Next()) {
        const auto& arc = ait.Value();
        auto it = M.find(arc.nextstate);
        if (it == M.end()) {
          M.emplace_hint(it, arc.nextstate, l + 1);
        } else {
          it->second = std::max(it->second, l + 1);
        }
      }
    }

    int max_length = 0;
    for (const auto& kv : M) {
      if (fst->Final(kv.first) != Weight::Zero()) {
        max_length = std::max(max_length, kv.second);
      }
    }
    return max_length;
  } else {
    return std::numeric_limits<int>::min();
  }
}

template<typename FST, typename F>
int ComputeMaxSubpathLength(FST* fst, const F& f, bool use_input = true) {
  typedef typename FST::Arc Arc;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  auto c_eps = f(0);
  typedef decltype(c_eps)  ClassType;

  if (fst->Properties(kAcyclic, true) & kAcyclic) {
    Connect(fst);
    TopSort(fst);
    if (fst->Start() == kNoStateId) { return -1; }

    std::map< StateId, std::map<ClassType, int> > M;
    M[fst->Start()][c_eps] = 0;
    for (StateIterator<FST> sit(*fst); !sit.Done(); sit.Next()) {
      const auto& s = sit.Value();
      const std::map<ClassType, int>& sm = M.find(s)->second;
      for (ArcIterator<FST> ait(*fst, s); !ait.Done(); ait.Next()) {
        const auto& arc = ait.Value();
        const auto& c_arc = f(use_input ? arc.ilabel : arc.olabel);
        auto it = M.find(arc.nextstate);
        if (it == M.end()) {
          it = M.emplace_hint(it, arc.nextstate, std::map<ClassType, int>());
        }
        std::map<ClassType, int>& sm2 = it->second;

        if (c_arc == c_eps) {
          for (const auto& kv : sm) {
            MapSetMax(sm2, kv.first, kv.second + 1);
          }
        } else {
          const auto max_prev = std::max(MapGetDefault(sm, c_arc, 0),
                                         MapGetDefault(sm, c_eps, 0));
          MapSetMax(sm2, c_arc, max_prev + 1);
        }
      }
    }

    int max_length = 0;
    for (const auto& kv1 : M) {
      for (const auto& kv2 : kv1.second) {
        max_length = std::max(max_length, kv2.second);
      }
    }
    return max_length;
  } else {
    return std::numeric_limits<int>::min();
  }
}

struct FstSummaryAcc {
  uint64_t num_fsts;
  uint64_t num_expanded;
  uint64_t num_mutable;
  uint64_t num_error;
  uint64_t num_acceptor;
  uint64_t num_idet;
  uint64_t num_odet;
  uint64_t num_isorted;
  uint64_t num_osorted;
  uint64_t num_weighted;
  uint64_t num_cyclic;
  uint64_t num_icyclic;
  uint64_t num_topsorted;
  double num_states;
  double num_arcs;
  double num_ioeps;
  double num_ieps;
  double num_oeps;
  double num_acc;
  double num_coacc;
  double num_cs;
  double num_cc;
  double num_scc;
  double num_ilm;
  double num_olm;
  double sum_sqr_states;
  double sum_sqr_arcs;
  double sum_sqr_ioeps;
  double sum_sqr_ieps;
  double sum_sqr_oeps;
  double sum_sqr_acc;
  double sum_sqr_coacc;
  double sum_sqr_cs;
  double sum_sqr_cc;
  double sum_sqr_scc;
  double sum_sqr_ilm;
  double sum_sqr_olm;

  long double num_paths;
  long double sum_sqr_paths;
  uint64_t num_inf_paths;

  int max_path_length;
  int max_subpath_length;

  FstSummaryAcc() :
      num_fsts(0),
      num_expanded(0),
      num_mutable(0),
      num_error(0),
      num_acceptor(0),
      num_idet(0),
      num_odet(0),
      num_isorted(0),
      num_osorted(0),
      num_weighted(0),
      num_cyclic(0),
      num_icyclic(0),
      num_topsorted(0),
      num_states(0),
      num_arcs(0),
      num_ioeps(0),
      num_ieps(0),
      num_oeps(0),
      num_acc(0),
      num_coacc(0),
      num_cs(0),
      num_cc(0),
      num_scc(0),
      num_ilm(0),
      num_olm(0),
      sum_sqr_states(0),
      sum_sqr_arcs(0),
      sum_sqr_ioeps(0),
      sum_sqr_ieps(0),
      sum_sqr_oeps(0),
      sum_sqr_acc(0),
      sum_sqr_coacc(0),
      sum_sqr_cs(0),
      sum_sqr_cc(0),
      sum_sqr_scc(0),
      sum_sqr_ilm(0),
      sum_sqr_olm(0),
      num_paths(0),
      sum_sqr_paths(0),
      num_inf_paths(0),
      max_path_length(std::numeric_limits<int>::min()),
      max_subpath_length(std::numeric_limits<int>::min()) {}
  friend std::ostream &operator<<(std::ostream &os,
                                  const FstSummaryAcc &summary) {
    const double N = summary.num_fsts;
    const double avg_states = N > 0 ? summary.num_states / N : 0;
    const double avg_arcs = N > 0 ? summary.num_arcs / N : 0;
    const double avg_ioeps = N > 0 ? summary.num_ioeps / N : 0;
    const double avg_ieps = N > 0 ? summary.num_ieps / N : 0;
    const double avg_oeps = N > 0 ? summary.num_oeps / N : 0;
    const double avg_acc = N > 0 ? summary.num_acc / N : 0;
    const double avg_coacc = N > 0 ? summary.num_coacc / N : 0;
    const double avg_cs = N > 0 ? summary.num_cs / N : 0;
    const double avg_cc = N > 0 ? summary.num_cc / N : 0;
    const double avg_scc = N > 0 ? summary.num_scc / N : 0;
    const double avg_ilm = N > 0 ? summary.num_ilm / N : 0;
    const double avg_olm = N > 0 ? summary.num_olm / N : 0;
    const long double avg_paths =
        (N - summary.num_inf_paths) > 0
        ? summary.num_paths / (N - summary.num_inf_paths) : 0;
    const double avg_expanded = N > 0 ? (100.0 * summary.num_expanded) / N : 0;
    const double avg_mutable = N > 0 ? (100.0 * summary.num_mutable) / N : 0;
    const double avg_error = N > 0 ? (100.0 * summary.num_error) / N : 0;
    const double avg_acceptor = N > 0 ? (100.0 * summary.num_acceptor) / N : 0;
    const double avg_idet = N > 0 ? (100.0 * summary.num_idet) / N : 0;
    const double avg_odet = N > 0 ? (100.0 * summary.num_odet) / N : 0;
    const double avg_isorted = N > 0 ? (100.0 * summary.num_isorted) / N : 0;
    const double avg_osorted = N > 0 ? (100.0 * summary.num_osorted) / N : 0;
    const double avg_weighted = N > 0 ? (100.0 * summary.num_weighted) / N : 0;
    const double avg_cyclic = N > 0 ? (100.0 * summary.num_cyclic) / N : 0;
    const double avg_icyclic = N > 0 ? (100.0 * summary.num_icyclic) / N : 0;
    const double avg_topsorted = N > 0 ? (100.0 * summary.num_topsorted) / N : 0;

    os << std::left;
    os << std::setw(50) << "# FSTs " << summary.num_fsts << std::endl;
    os << std::setw(50) << "avg. of states" << avg_states << std::endl;
    os << std::setw(50) << "avg. of arcs" << avg_arcs << std::endl;
    os << std::setw(50) << "avg. of input/output epsilons" << avg_ioeps
       << std::endl;
    os << std::setw(50) << "avg. of input epsilons" << avg_ieps << std::endl;
    os << std::setw(50) << "avg. of output epsilons" << avg_oeps << std::endl;
    os << std::setw(50) << "avg. of accessible states" << avg_acc << std::endl;
    os << std::setw(50) << "avg. of coaccessible states" << avg_coacc
       << std::endl;
    os << std::setw(50) << "avg. of connected states" << avg_cs << std::endl;
    os << std::setw(50) << "avg. of connected components" << avg_cc
       << std::endl;
    os << std::setw(50) << "avg. of strongly conn components" << avg_scc
       << std::endl;
    os << std::setw(50) << "avg. of paths" << avg_paths << std::endl;
    os << std::setw(50) << "avg. input label multiplicity" << avg_ilm
       << std::endl;
    os << std::setw(50) << "avg. output label multiplicity" << avg_olm
       << std::endl;
    if (summary.max_path_length >= 0) {
      std::cout << std::setw(50) << "max. path length"
                << summary.max_path_length << std::endl;
    } else {
      std::cout << std::setw(50) << "max. path length"
                << "none" << std::endl;
    }
    if (summary.max_path_length >= 0) {
      std::cout << std::setw(50) << "max. subpath length"
                << summary.max_subpath_length << std::endl;
    } else {
      std::cout << std::setw(50) << "max. subpath length"
                << "none" << std::endl;
    }
    os << std::setw(50) << "% expanded" << avg_expanded << std::endl;
    os << std::setw(50) << "% mutable" << avg_mutable << std::endl;
    os << std::setw(50) << "% error" << avg_error << std::endl;
    os << std::setw(50) << "% acceptor" << avg_acceptor << std::endl;
    os << std::setw(50) << "% input deterministic" << avg_idet << std::endl;
    os << std::setw(50) << "% output deterministic" << avg_odet << std::endl;
    os << std::setw(50) << "% input label sorted" << avg_isorted << std::endl;
    os << std::setw(50) << "% output label sorted" << avg_osorted << std::endl;
    os << std::setw(50) << "% weighted" << avg_weighted << std::endl;
    os << std::setw(50) << "% cyclic" << avg_cyclic << std::endl;
    os << std::setw(50) << "% cyclic at initial state" << avg_icyclic
       << std::endl;
    os << std::setw(50) << "% top sorted" << avg_topsorted << std::endl;
    return os;
  }
};

template<typename FstReader, typename F>
void UpdateFstSummary(
    const std::string &rspecifier, FstSummaryAcc *acc,
    const F* func = nullptr) {
  typedef typename FstReader::T FST;
  for (FstReader fst_reader(rspecifier); !fst_reader.Done();
       fst_reader.Next()) {
    FST f = fst_reader.Value();
    fst_reader.FreeCurrent();
    FstInfo fst_info(f, true);

    acc->num_fsts++;
    acc->num_states += fst_info.NumStates();
    acc->num_arcs += fst_info.NumArcs();
    acc->num_ioeps += fst_info.NumEpsilons();
    acc->num_ieps += fst_info.NumInputEpsilons();
    acc->num_oeps += fst_info.NumOutputEpsilons();
    acc->num_acc += fst_info.NumAccessible();
    acc->num_coacc += fst_info.NumCoAccessible();
    acc->num_cs += fst_info.NumConnected();
    acc->num_cc += fst_info.NumCc();
    acc->num_scc += fst_info.NumScc();
    acc->num_ilm += fst_info.InputLabelMultiplicity();
    acc->num_olm += fst_info.OutputLabelMultiplicity();
    acc->sum_sqr_states += fst_info.NumStates() * fst_info.NumStates();
    acc->sum_sqr_arcs += fst_info.NumArcs() * fst_info.NumArcs();
    acc->sum_sqr_ioeps += fst_info.NumEpsilons() * fst_info.NumEpsilons();
    acc->sum_sqr_ieps +=
        fst_info.NumInputEpsilons() * fst_info.NumInputEpsilons();
    acc->sum_sqr_oeps +=
        fst_info.NumOutputEpsilons() * fst_info.NumOutputEpsilons();
    acc->sum_sqr_acc += fst_info.NumAccessible() * fst_info.NumAccessible();
    acc->sum_sqr_coacc +=
        fst_info.NumCoAccessible() * fst_info.NumCoAccessible();
    acc->sum_sqr_cs += fst_info.NumConnected() * fst_info.NumConnected();
    acc->sum_sqr_cc += fst_info.NumCc() * fst_info.NumCc();
    acc->sum_sqr_scc += fst_info.NumScc() * fst_info.NumScc();
    acc->sum_sqr_ilm +=
        fst_info.InputLabelMultiplicity() * fst_info.InputLabelMultiplicity();
    acc->sum_sqr_olm +=
        fst_info.OutputLabelMultiplicity() * fst_info.OutputLabelMultiplicity();

    const long double num_paths = ComputeNumberOfPaths<FST>(&f);
    if (num_paths < std::numeric_limits<long double>::infinity()) {
      acc->num_paths += num_paths;
      acc->sum_sqr_paths += num_paths * num_paths;
    } else {
      acc->num_inf_paths++;
    }

    // Print info about the max length of all paths
    acc->max_path_length = std::max(acc->max_path_length,
                                    ComputeMaxPathLength<FST>(&f));
    if (func != nullptr) {
      acc->max_subpath_length = std::max(
          acc->max_subpath_length,
          ComputeMaxSubpathLength<FST, F>(&f, *func));
    }

    if (fst_info.Properties() & fst::kExpanded)
      acc->num_expanded++;
    if (fst_info.Properties() & fst::kMutable)
      acc->num_mutable++;
    if (fst_info.Properties() & fst::kError)
      acc->num_error++;
    if (fst_info.Properties() & fst::kAcceptor)
      acc->num_acceptor++;
    if (fst_info.Properties() & fst::kIDeterministic)
      acc->num_idet++;
    if (fst_info.Properties() & fst::kODeterministic)
      acc->num_odet++;
    if (fst_info.Properties() & fst::kILabelSorted)
      acc->num_isorted++;
    if (fst_info.Properties() & fst::kOLabelSorted)
      acc->num_osorted++;
    if (fst_info.Properties() & fst::kWeighted)
      acc->num_weighted++;
    if (fst_info.Properties() & fst::kCyclic)
      acc->num_cyclic++;
    if (fst_info.Properties() & fst::kInitialCyclic)
      acc->num_icyclic++;
    if (fst_info.Properties() & fst::kTopSorted)
      acc->num_topsorted++;
  }
}

template<typename FstReader, typename F>
void PrintFstInfo(const std::string &rspecifier, const F* func) {
  typedef typename FstReader::T FST;
  typedef typename FST::Arc Arc;
  for (FstReader fst_reader(rspecifier); !fst_reader.Done();
       fst_reader.Next()) {
    const std::string key = fst_reader.Key();
    FST f = fst_reader.Value();
    fst_reader.FreeCurrent();
    fst::FstInfo fst_info(f, true);
    std::cout << std::left << key << std::endl;
    std::cout << std::setw(50) << "# of states" << fst_info.NumStates()
              << std::endl;
    std::cout << std::setw(50) << "# of arcs" << fst_info.NumArcs()
              << std::endl;
    std::cout << std::setw(50) << "initial state" << fst_info.Start()
              << std::endl;
    std::cout << std::setw(50) << "# of input/output epsilons"
              << fst_info.NumEpsilons() << std::endl;
    std::cout << std::setw(50) << "# of input epsilons"
              << fst_info.NumInputEpsilons() << std::endl;
    std::cout << std::setw(50) << "# of output epsilons"
              << fst_info.NumOutputEpsilons() << std::endl;
    std::cout << std::setw(50) << "# of accessible states"
              << fst_info.NumAccessible() << std::endl;
    std::cout << std::setw(50) << "# of coaccessible states"
              << fst_info.NumCoAccessible() << std::endl;
    std::cout << std::setw(50) << "# of connected states"
              << fst_info.NumConnected() << std::endl;
    std::cout << std::setw(50) << "# of connected components"
              << fst_info.NumCc() << std::endl;
    std::cout << std::setw(50) << "# of strongly conn components"
              << fst_info.NumScc() << std::endl;
    // Print number of paths
    const long double num_paths = ComputeNumberOfPaths<FST>(&f);
    if (num_paths > std::pow(2, std::numeric_limits<long double>::digits)) {
      std::cout << std::setw(50) << "# of paths" << num_paths << std::endl;
    } else {
      std::cout << std::setw(50) << "# of paths"
                << static_cast<uint64_t>(num_paths) << std::endl;
    }
    std::cout << std::setw(50) << "input label multiplicity"
              << fst_info.InputLabelMultiplicity() << std::endl;
    std::cout << std::setw(50) << "output label multiplicity"
              << fst_info.OutputLabelMultiplicity() << std::endl;
    // Print info about the max length of all paths
    const int max_path_length = ComputeMaxPathLength<FST>(&f);
    if (max_path_length >= 0) {
      std::cout << std::setw(50) << "max. path length"
                << max_path_length << std::endl;
    } else {
      // Fst is empty or not acyclic
      std::cout << std::setw(50) << "max. path length"
                << "none" << std::endl;
    }
    if (func != nullptr) {
      const int max_subpath_length = ComputeMaxSubpathLength<FST, F>(&f, *func);
      if (max_subpath_length >= 0) {
        std::cout << std::setw(50) << "max. subpath length"
                  << max_subpath_length << std::endl;
      } else {
        // Fst is empty or not acyclic
        std::cout << std::setw(50) << "max. subpath length"
                  << "none" << std::endl;
      }
    } else {
      std::cout << std::setw(50) << "max. subpath length"
                << "none" << std::endl;
    }

    uint64 prop = 1;
    for (int i = 0; i < 64; ++i, prop <<= 1) {
      if (prop & fst::kBinaryProperties) {
        char value = 'n';
        if (fst_info.Properties() & prop) value = 'y';
        std::cout << std::setw(50) << fst::PropertyNames[i] << value
                  << std::endl;
      } else if (prop & fst::kPosTrinaryProperties) {
        char value = '?';
        if (fst_info.Properties() & prop) value = 'y';
        else if (fst_info.Properties() & prop << 1) value = 'n';
        std::cout << std::setw(50) << fst::PropertyNames[i] << value
                  << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

}   // namespace fst

#endif //KALDI_LATTICE_UTILS_FST_INFO_H_

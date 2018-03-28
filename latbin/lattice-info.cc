// latbin/lattice-info.cc

// Copyright 2009-2011  Microsoft Corporation
//                2013  Johns Hopkins University (author: Daniel Povey)

// See ../../COPYING for clarification regarding multiple authors
//
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


#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/fstext-lib.h"
#include "lat/kaldi-lattice.h"
#include "fst/script/info-impl.h"

#include <iomanip>

namespace kaldi {

template <typename Lattice>
long double ComputeNumberOfPaths(Lattice* lat) {
  const long double INF = std::numeric_limits<long double>::infinity();
  if (lat->Properties(fst::kAcyclic, true) & fst::kAcyclic) {
    fst::TopSort(lat);
    // (Log)number of paths arriving to each state.
    std::vector<long double> num_paths(lat->NumStates(), 0);
    num_paths[lat->Start()] = 1;
    // (Log)number of full paths.
    long double num_full_paths = 0;
    for (fst::StateIterator<Lattice> siter(*lat); !siter.Done(); siter.Next()) {
      const typename Lattice::StateId s = siter.Value();
      for (fst::ArcIterator<Lattice> aiter(*lat, s); !aiter.Done();
           aiter.Next()) {
        const typename Lattice::Arc arc = aiter.Value();
        num_paths[arc.nextstate] += num_paths[s];
      }
      if (lat->Final(s) != Lattice::Arc::Weight::Zero()) {
        num_full_paths += num_paths[s];
      }
    }
    return num_full_paths;
  } else {
    return INF;
  }
}

struct LatticeSummaryAcc {
  uint64_t num_lattices;
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

  LatticeSummaryAcc() :
      num_lattices(0),
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
      num_inf_paths(0)
  { }
  friend std::ostream& operator<<(std::ostream& os,
                                  const LatticeSummaryAcc& summary) {
    const double avg_states = summary.num_lattices > 0 ?
        summary.num_states / summary.num_lattices : 0;
    const double avg_arcs = summary.num_lattices > 0 ?
        summary.num_arcs / summary.num_lattices : 0;
    const double avg_ioeps = summary.num_lattices > 0 ?
        summary.num_ioeps / summary.num_lattices : 0;
    const double avg_ieps = summary.num_lattices > 0 ?
        summary.num_ieps / summary.num_lattices : 0;
    const double avg_oeps = summary.num_lattices > 0 ?
        summary.num_oeps / summary.num_lattices : 0;
    const double avg_acc = summary.num_lattices > 0 ?
        summary.num_acc / summary.num_lattices : 0;
    const double avg_coacc = summary.num_lattices > 0 ?
        summary.num_coacc / summary.num_lattices : 0;
    const double avg_cs = summary.num_lattices > 0 ?
        summary.num_cs / summary.num_lattices : 0;
    const double avg_cc = summary.num_lattices > 0 ?
        summary.num_cc / summary.num_lattices : 0;
    const double avg_scc = summary.num_lattices > 0 ?
        summary.num_scc / summary.num_lattices : 0;
    const double avg_ilm = summary.num_lattices > 0 ?
        summary.num_ilm / summary.num_lattices : 0;
    const double avg_olm = summary.num_lattices > 0 ?
        summary.num_olm / summary.num_lattices : 0;
    const long double avg_paths = summary.num_lattices > 0 ?
        summary.num_paths / summary.num_lattices : 0;
    const double avg_expanded = summary.num_lattices > 0 ?
        100.0 * summary.num_expanded / (double)summary.num_lattices : 0;
    const double avg_mutable = summary.num_lattices > 0 ?
        100.0 * summary.num_mutable / (double)summary.num_lattices : 0;
    const double avg_error = summary.num_lattices > 0 ?
        100.0 * summary.num_error / (double)summary.num_lattices : 0;
    const double avg_acceptor = summary.num_lattices > 0 ?
        100.0 * summary.num_acceptor / (double)summary.num_lattices : 0;
    const double avg_idet = summary.num_lattices > 0 ?
        100.0 * summary.num_idet / (double)summary.num_lattices : 0;
    const double avg_odet = summary.num_lattices > 0 ?
        100.0 * summary.num_odet / (double)summary.num_lattices : 0;
    const double avg_isorted = summary.num_lattices > 0 ?
        100.0 * summary.num_isorted / (double)summary.num_lattices : 0;
    const double avg_osorted = summary.num_lattices > 0 ?
        100.0 * summary.num_osorted / (double)summary.num_lattices : 0;
    const double avg_weighted = summary.num_lattices > 0 ?
        100.0 * summary.num_weighted / (double)summary.num_lattices : 0;
    const double avg_cyclic = summary.num_lattices > 0 ?
        100.0 * summary.num_cyclic / (double)summary.num_lattices : 0;
    const double avg_icyclic = summary.num_lattices > 0 ?
        100.0 * summary.num_icyclic / (double)summary.num_lattices : 0;
    const double avg_topsorted = summary.num_lattices > 0 ?
        100.0 * summary.num_topsorted / (double)summary.num_lattices : 0;
    os << std::left;
    os << std::setw(50) << "# lattices " << summary.num_lattices << std::endl;
    os << std::setw(50) << "avg. of states" << avg_states << std::endl;
    os << std::setw(50) << "avg. of arcs" << avg_arcs << std::endl;
    os << std::setw(50) << "avg. of input/output epsilons" << avg_ioeps << std::endl;
    os << std::setw(50) << "avg. of input epsilons" << avg_ieps << std::endl;
    os << std::setw(50) << "avg. of output epsilons" << avg_oeps << std::endl;
    os << std::setw(50) << "avg. of accessible states" << avg_acc << std::endl;
    os << std::setw(50) << "avg. of coaccessible states" << avg_coacc << std::endl;
    os << std::setw(50) << "avg. of connected states" << avg_cs << std::endl;
    os << std::setw(50) << "avg. of connected components" << avg_cc << std::endl;
    os << std::setw(50) << "avg. of strongly conn components" << avg_scc << std::endl;
    os << std::setw(50) << "avg. of paths" << avg_paths << std::endl;
    os << std::setw(50) << "% input label multiplicity" << avg_ilm << std::endl;
    os << std::setw(50) << "% output label multiplicity" << avg_olm << std::endl;
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
    os << std::setw(50) << "% cyclic at initial state" << avg_icyclic << std::endl;
    os << std::setw(50) << "% top sorted" << avg_topsorted << std::endl;
    return os;
  }
};

template <typename LatticeReader>
void UpdateLatticeSummary(const std::string& lattice_rspecifier, LatticeSummaryAcc* acc) {
  typedef typename LatticeReader::T Lattice;
  for (LatticeReader lattice_reader(lattice_rspecifier);
       !lattice_reader.Done(); lattice_reader.Next()) {
    Lattice lat = lattice_reader.Value();
    lattice_reader.FreeCurrent();
    fst::FstInfo fst_info(lat, true);
    const long double num_paths = ComputeNumberOfPaths<Lattice>(&lat);
    acc->num_lattices++;
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
    acc->sum_sqr_ieps += fst_info.NumInputEpsilons() * fst_info.NumInputEpsilons();
    acc->sum_sqr_oeps += fst_info.NumOutputEpsilons() * fst_info.NumOutputEpsilons();
    acc->sum_sqr_acc += fst_info.NumAccessible() * fst_info.NumAccessible();
    acc->sum_sqr_coacc += fst_info.NumCoAccessible() * fst_info.NumCoAccessible();
    acc->sum_sqr_cs += fst_info.NumConnected() * fst_info.NumConnected();
    acc->sum_sqr_cc += fst_info.NumCc() * fst_info.NumCc();
    acc->sum_sqr_scc += fst_info.NumScc() * fst_info.NumScc();
    acc->sum_sqr_ilm += fst_info.InputLabelMultiplicity() * fst_info.InputLabelMultiplicity();
    acc->sum_sqr_olm += fst_info.OutputLabelMultiplicity() * fst_info.OutputLabelMultiplicity();

    if (num_paths < std::numeric_limits<long double>::infinity()) {
      acc->num_paths += num_paths;
      acc->sum_sqr_paths += num_paths * num_paths;
    } else {
      acc->num_inf_paths++;
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
      acc->num_isorted++;
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

template <typename LatticeReader>
void PrintLatticeInfo(const std::string& lattice_rspecifier) {
  typedef typename LatticeReader::T Lattice;
  for (LatticeReader lattice_reader(lattice_rspecifier);
       !lattice_reader.Done(); lattice_reader.Next()) {
    const std::string key = lattice_reader.Key();
    Lattice lat = lattice_reader.Value();
    lattice_reader.FreeCurrent();
    fst::FstInfo fst_info(lat, true);
    const long double num_paths = ComputeNumberOfPaths<Lattice>(&lat);
    std::cout << std::left << key << std::endl;
    std::cout << std::setw(50) << "# of states" << fst_info.NumStates() << std::endl;
    std::cout << std::setw(50) << "# of arcs" << fst_info.NumArcs() << std::endl;
    std::cout << std::setw(50) << "initial state" << fst_info.Start() << std::endl;
    std::cout << std::setw(50) << "# of input/output epsilons" << fst_info.NumEpsilons() << std::endl;
    std::cout << std::setw(50) << "# of input epsilons" << fst_info.NumInputEpsilons() << std::endl;
    std::cout << std::setw(50) << "# of output epsilons" << fst_info.NumOutputEpsilons() << std::endl;
    std::cout << std::setw(50) << "# of accessible states" << fst_info.NumAccessible() << std::endl;
    std::cout << std::setw(50) << "# of coaccessible states" << fst_info.NumCoAccessible() << std::endl;
    std::cout << std::setw(50) << "# of connected states" << fst_info.NumConnected() << std::endl;
    std::cout << std::setw(50) << "# of connected components" << fst_info.NumCc() << std::endl;
    std::cout << std::setw(50) << "# of strongly conn components" << fst_info.NumScc() << std::endl;
    std::cout << std::setw(50) << "# of paths" << num_paths << std::endl;
    std::cout << std::setw(50) << "input label multiplicity" << fst_info.InputLabelMultiplicity() << std::endl;
    std::cout << std::setw(50) << "output label multiplicity" << fst_info.OutputLabelMultiplicity() << std::endl;
    uint64 prop = 1;
    for (int i = 0; i < 64; ++i, prop <<= 1) {
      if (prop & fst::kBinaryProperties) {
        char value = 'n';
        if (fst_info.Properties() & prop) value = 'y';
        std::cout << std::setw(50) << fst::PropertyNames[i] << value << std::endl;
      } else if (prop & fst::kPosTrinaryProperties) {
        char value = '?';
        if (fst_info.Properties() & prop) value = 'y';
        else if (fst_info.Properties() & prop << 1) value = 'n';
        std::cout << std::setw(50) << fst::PropertyNames[i] << value << std::endl;
      }
    }
    std::cout << std::endl;
  }
}


}  // namespace kaldi

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;

    const char *usage =
        "Prints out a summary of the lattices in a set of archives or\n"
        "the individual information about each of the lattices\n"
        "(--summary=false), similar to what fstinfo does.\n"
        "Usage: lattice-info [options] lattice-rspecifier1 [lattice-rspecifier2 ...]\n"
        " e.g.: lattice-info --summary=false ark:1.lats ark:2.lats\n";

    ParseOptions po(usage);
    bool compact = true;
    bool summary = true;
    std::string include_rxfilename;
    std::string exclude_rxfilename;

    po.Register("compact", &compact, "If true, work with lattices in compact form.");
    po.Register("summary", &summary, "If true, summarizes the information of all lattices.");
    po.Register("include", &include_rxfilename,
                "Text file, the first field of each "
                "line being interpreted as the "
                "utterance-id whose lattices will be included");
    po.Register("exclude", &exclude_rxfilename,
                "Text file, the first field of each "
                "line being interpreted as an utterance-id "
                "whose lattices will be excluded");

    po.Read(argc, argv);

    if (po.NumArgs() < 1) {
      po.PrintUsage();
      exit(1);
    }

    if (include_rxfilename != "" && exclude_rxfilename != "") {
      KALDI_ERR << "should not have both --exclude and --include option!";
    }

    LatticeSummaryAcc summary_acc;
    for (int32 a = 1; a <= po.NumArgs(); ++a) {
      const std::string lats_rspecifier = po.GetArg(a);
      if (compact) {
        if (summary) {
          UpdateLatticeSummary<SequentialCompactLatticeReader>(lats_rspecifier,
                                                               &summary_acc);
        } else {
          PrintLatticeInfo<SequentialCompactLatticeReader>(lats_rspecifier);
        }
      } else {
        if (summary) {
          UpdateLatticeSummary<SequentialLatticeReader>(lats_rspecifier,
                                                        &summary_acc);
        } else {
          PrintLatticeInfo<SequentialLatticeReader>(lats_rspecifier);
        }
      }
    }
    if (summary) {
      std::cout << summary_acc;
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

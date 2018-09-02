// latbin/lattice-prune-arcs.cc

//                2017  Joan Puigcerver

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
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

namespace kaldi {

template <typename T>
struct FirstElementTupleSorter {
  bool operator()(const T& a, const T& b) const {
    return std::get<0>(a) < std::get<0>(b);
  }
};

template<class LatType>
// could be Lattice or CompactLattice
void PruneLatticeArcs(double beam, LatType *lat) {
  typedef typename LatType::Arc Arc;
  typedef typename Arc::StateId StateId;

  if(lat->Start() == fst::kNoStateId) { return; }

  std::vector<double> alphas, betas;
  const double total = ComputeLatticeAlphasAndBetas(*lat, false, &alphas, &betas);
  const double cost_cutoff = beam - total;
  for (size_t i = 0; i < alphas.size(); ++i) {
    alphas[i] = -alphas[i];
    betas[i] = -betas[i];
  }

  std::vector< std::tuple<double, StateId, Arc> > arcs;
  for (fst::StateIterator<LatType> sit(*lat); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    for (fst::ArcIterator<LatType> ait(*lat, s); !ait.Done(); ait.Next()) {
        const auto& arc = ait.Value();
        const double cost_arc = ConvertToCost(arc.weight);
        const double cost_through = cost_arc + alphas[s] + betas[arc.nextstate];
        arcs.emplace_back(cost_through, s, arc);
    }
    lat->DeleteArcs(s);
  }
  // Sort arcs in increasing order of probability.
  FirstElementTupleSorter< std::tuple<double, StateId, Arc> > sorter;
  std::sort(arcs.begin(), arcs.end(), sorter);


  size_t i = 0;
  double cost_acc = std::numeric_limits<double>::infinity();
  for (; i < arcs.size(); ++i) {
    const auto& tup = arcs[i];
    cost_acc = -LogAdd(-cost_acc, -std::get<0>(tup));
    if (cost_acc < cost_cutoff) break;
  }

  if (i == arcs.size()) {
    lat->DeleteStates();
  }

  for (; i < arcs.size(); ++i) {
    const auto& tup = arcs[i];
    lat->AddArc(std::get<1>(tup), std::get<2>(tup));
  }

  fst::Connect(lat);
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
        "Iteratively reduce the beam of the lattice until a maximum number "
            "of arcs and states is achieved.";

    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    ParseOptions po(usage);
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("beam", &beam, "");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Read(argc, argv);

    if (po.NumArgs() < 2) {
      po.PrintUsage();
      exit(1);
    }

    if (beam <= 0.0) {
      KALDI_ERR << "--beam_ratio must be in the open range (0.0, inf).";
    }

    const auto lattice_scale = fst::LatticeScale(graph_scale, acoustic_scale);
    const auto inv_lattice_scale =
        fst::LatticeScale(1.0 / graph_scale, 1.0 / acoustic_scale);
    SequentialCompactLatticeReader lattice_reader(po.GetArg(1));
    CompactLatticeWriter lattice_writer(po.GetArg(2));
    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      // Read input lattice
      const std::string lattice_key = lattice_reader.Key();
      CompactLattice &ilat = lattice_reader.Value();

      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(lattice_scale, &ilat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(insertion_penalty, &ilat);
      // Topological order of nodes in the lattice
      TopSortCompactLatticeIfNeeded(&ilat);

      const int32 original_num_arcs = NumArcs(ilat);
      const int32 original_num_states = ilat.NumStates();
      PruneLatticeArcs(beam, &ilat);
      const int32 num_arcs = NumArcs(ilat);
      const int32 num_states = ilat.NumStates();

      // Lattice is written in the original scale
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(inv_lattice_scale, &ilat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(-insertion_penalty, &ilat);

      lattice_writer.Write(lattice_key, ilat);
      KALDI_LOG << "Lattice " << lattice_key << " pruned #states from "
                << original_num_states << " to " << num_states
                << " and #arcs from " << original_num_arcs << " to "
                << num_arcs;
    }

    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

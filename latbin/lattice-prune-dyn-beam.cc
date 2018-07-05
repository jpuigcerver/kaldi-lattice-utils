// latbin/lattice-prune-dyn-beam.cc

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

template<class LatType>
// could be Lattice or CompactLattice
double ComputeLatticeBeam(const LatType &lat) {
  typedef typename LatType::Arc Arc;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId StateId;
  const int32 num_states = lat.NumStates();
  if (num_states == 0) return 0.0;
  std::vector<double> forward_cost(lat.NumStates(),
                                   std::numeric_limits<double>::infinity());
  // Note: If FST is topologically ordered, then state 0 is the start
  forward_cost[lat.Start()] = 0.0;
  // Compute Viterbi forward
  double best_final_cost = std::numeric_limits<double>::infinity();
  for (int32 state = 0; state < num_states; state++) {
    double this_forward_cost = forward_cost[state];
    for (fst::ArcIterator <LatType> aiter(lat, state);
         !aiter.Done();
         aiter.Next()) {
      const Arc &arc(aiter.Value());
      StateId nextstate = arc.nextstate;
      KALDI_ASSERT(nextstate > state && nextstate < num_states);
      double next_forward_cost = this_forward_cost +
          ConvertToCost(arc.weight);
      if (forward_cost[nextstate] > next_forward_cost)
        forward_cost[nextstate] = next_forward_cost;
    }
    Weight final_weight = lat.Final(state);
    double this_final_cost = this_forward_cost +
        ConvertToCost(final_weight);
    if (this_final_cost < best_final_cost)
      best_final_cost = this_final_cost;
  }
  double cutoff = best_final_cost;

  // Go backwards updating the backward probs (which share memory with the
  // forward probs), and updating the cutoff;
  std::vector<double> &backward_cost(forward_cost);
  for (int32 state = num_states - 1; state >= 0; state--) {
    double this_forward_cost = forward_cost[state];
    double this_backward_cost = ConvertToCost(lat.Final(state));
    if (this_backward_cost + this_forward_cost > cutoff
        && this_backward_cost != std::numeric_limits<double>::infinity())
      cutoff = this_backward_cost + this_forward_cost;
    for (fst::ArcIterator <LatType> aiter(lat, state);
         !aiter.Done();
         aiter.Next()) {
      Arc arc(aiter.Value());
      StateId nextstate = arc.nextstate;
      KALDI_ASSERT(nextstate > state && nextstate < num_states);
      double arc_cost = ConvertToCost(arc.weight),
          arc_backward_cost = arc_cost + backward_cost[nextstate],
          this_fb_cost = this_forward_cost + arc_backward_cost;
      if (arc_backward_cost < this_backward_cost)
        this_backward_cost = arc_backward_cost;
      if (this_fb_cost > cutoff) {
        cutoff = this_fb_cost;
      }
    }
    backward_cost[state] = this_backward_cost;
  }
  // beam = cutoff - best_final_cost
  return cutoff - best_final_cost;
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

    int32 max_num_arcs = std::numeric_limits<int32>::max();
    int32 max_num_states = std::numeric_limits<int32>::max();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam_ratio = 0.9;
    BaseFloat min_beam = 1e-3;
    ParseOptions po(usage);
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("beam-ratio", &beam_ratio,
                "Reduce the maximum beam by this ratio at each iteration.");
    po.Register("min-beam", &min_beam,
                "Minimum beam threshold");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                    "label (typically, equivalent to word insertion penalty).");
    po.Register("max-arcs", &max_num_arcs,
                "Maximum number of arcs of each lattice.");
    po.Register("max-states", &max_num_states,
                "Maximum number of states of each lattice.");
    po.Read(argc, argv);

    if (po.NumArgs() < 2) {
      po.PrintUsage();
      exit(1);
    }

    if (beam_ratio <= 0.0 || beam_ratio >= 1.0) {
      KALDI_ERR << "--beam_ratio must be in the open range (0.0, 1.0).";
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

      const auto original_beam = ComputeLatticeBeam(ilat);
      const int32 original_num_arcs = NumArcs(ilat);
      const int32 original_num_states = ilat.NumStates();

      auto beam = original_beam;
      int32 num_arcs = original_num_arcs;
      int32 num_states = original_num_states;
      for (size_t n_try = 1; beam > min_beam &&
          (num_arcs > max_num_arcs || num_states > max_num_states);
           ++n_try) {
        beam = beam_ratio * beam;
        PruneLattice(beam, &ilat);
        const auto t_num_arcs = NumArcs(ilat);
        const auto t_num_states = ilat.NumStates();
        KALDI_VLOG(1) << "Lattice " << lattice_key << " pruned with beam = "
                      << beam << " (" << n_try << " trial): "
                      << " pruned #states from " << num_states
                      << " to " << t_num_states
                      << " and #arcs from " << num_arcs
                      << " to " << t_num_arcs;
        num_arcs = t_num_arcs;
        num_states = t_num_states;
      }

      // Lattice is written in the original scale
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(inv_lattice_scale, &ilat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(-insertion_penalty, &ilat);

      lattice_writer.Write(lattice_key, ilat);
      if (original_num_states == num_states &&
          original_num_arcs == num_arcs) {
        KALDI_LOG << "Lattice " << lattice_key << " was not pruned (beam = "
                  << original_beam << ", # states = " << original_num_states
                  << ", # arcs = " << original_num_arcs << ")";
      } else {
        KALDI_LOG << "Lattice " << lattice_key << " pruned #states from "
                  << original_num_states << " to " << num_states
                  << " and #arcs from " << original_num_arcs << " to "
                  << num_arcs << " (beam reduced from " << original_beam
                  << " to " << beam << ")";
      }
    }

    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

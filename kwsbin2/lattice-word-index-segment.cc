// kwsbin2/lattice-word-index-segment.cc

// Copyright (c) 2016-2018 Joan Puigcerver <joapuipe@upv.es>

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

#include <algorithm>
#include <string>
#include <tuple>

#include "base/kaldi-common.h"
#include "base/timer.h"
#include "fstext/kaldi-fst-io.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"


namespace kaldi {
void ProcessLattice(
    const std::string &key, CompactLattice *clat, BaseFloat acoustic_scale,
    BaseFloat graph_scale, BaseFloat insertion_penalty, BaseFloat beam,
    std::vector<int32> *state_times,
    std::vector<double> *fw, std::vector<double> *bw, double *total_lkh) {
  const int64 narcs = NumArcs(*clat);
  const int64 nstates = clat->NumStates();
  // Acoustic scale
  if (acoustic_scale != 1.0 || graph_scale != 1.0) {
    ScaleLattice(fst::LatticeScale(graph_scale, acoustic_scale), clat);
  }
  // Word insertion penalty
  if (insertion_penalty != 0.0) {
    AddWordInsPenToCompactLattice(insertion_penalty, clat);
  }
  // Lattice pruning
  if (beam != std::numeric_limits<BaseFloat>::infinity()) {
    PruneLattice(beam, clat);
  }
  const int64 pruned_narcs = NumArcs(*clat);
  const int64 pruned_nstates = clat->NumStates();
  KALDI_VLOG(1) << "Lattice " << key << ": pruned #states from "
                << nstates << " to " << pruned_nstates << " and #arcs from "
                << narcs << " to " << pruned_narcs;

  if (clat->Start() != fst::kNoStateId) {
    // If needed, sort the compact lattice in topological order
    TopSortCompactLatticeIfNeeded(clat);
    CompactLatticeStateTimes(*clat, state_times);
    *total_lkh = ComputeLatticeAlphasAndBetas(*clat, false, fw, bw);
  } else {
    state_times->clear();
    fw->clear();
    bw->clear();
    *total_lkh = 0;
  }
}

class LatticeScorerTask {
 public:
  typedef std::map<int32, std::map<std::tuple<int32, int32>, double>> MMT;
  typedef TableWriter<BasicTupleVectorHolder<int32, int32, int32, double>>
      SegmentScoreWriter;

  LatticeScorerTask(
      const std::string &key, const CompactLattice &clat,
      const std::set<int32> &exclude_labels,
      BaseFloat acoustic_scale, BaseFloat graph_scale,
      BaseFloat insertion_penalty, BaseFloat beam,
      SegmentScoreWriter* score_writer) :
      key_(key), clat_(clat), exclude_labels_(exclude_labels),
      acoustic_scale_(acoustic_scale), graph_scale_(graph_scale),
      insertion_penalty_(insertion_penalty), beam_(beam),
      score_writer_(score_writer) {}

  ~LatticeScorerTask() {
    // Put all segments of all words into a single vector.
    std::vector<std::tuple<int32, int32, int32, double>> segment_scores;
    for (const auto& ws : accumulator_) {
      const int32 word_label = ws.first;
      const auto& segments = ws.second;
      for (const auto& ttp : segments) {
        const auto& tt = ttp.first;
        const double lkh = ttp.second;
        segment_scores.emplace_back(
            word_label, std::get<0>(tt), std::get<1>(tt), lkh - total_lkh_);
      }
    }
    // Sort the vector according to:
    //   1) Descending probability
    //   2) Ascending word label
    //   3) Ascending t0
    //   4) Ascending t1
    std::sort(segment_scores.begin(), segment_scores.end(),
              [](const std::tuple<int32, int32, int32, double> &a,
                 const std::tuple<int32, int32, int32, double> &b) -> bool {
                if (std::get<3>(b) != std::get<3>(a)) {
                  return std::get<3>(b) < std::get<3>(a);
                } else if (std::get<0>(a) != std::get<0>(b)) {
                  return std::get<0>(a) < std::get<0>(b);
                } else if (std::get<1>(a) != std::get<1>(b)) {
                  return std::get<1>(a) < std::get<1>(b);
                } else {
                  return std::get<2>(a) < std::get<2>(b);
                }
              });

    score_writer_->Write(key_, segment_scores);

    KALDI_VLOG(1) << "Lattice " << key_ << ": done in " << elapsed_
                  << " seconds.";
  }

  void operator()() {
    Timer timer;

    // Prepare lattice to compute the scores.
    ProcessLattice(
        key_, &clat_, acoustic_scale_, graph_scale_, insertion_penalty_,
        beam_, &state_times_, &fw_, &bw_, &total_lkh_);

    for (fst::StateIterator<CompactLattice> sit(clat_); !sit.Done();
         sit.Next()) {
      for (fst::ArcIterator<CompactLattice> ait(clat_, sit.Value());
           !ait.Done(); ait.Next()) {
        const auto &arc = ait.Value();
        const auto &label = arc.olabel;
        if (label != 0 && exclude_labels_.count(label) == 0) {
          // Initial and final timesteps of the arc.
          const int32 t0 = state_times_[sit.Value()];
          const int32 t1 = state_times_[arc.nextstate];
          // Likelhood of the arc.
          const double arc_lkh = -ConvertToCost(arc.weight);
          // Likelihood of all paths through the arc.
          const double lkh_through_arc =
              fw_[sit.Value()] + arc_lkh + bw_[arc.nextstate];
          // Add (or increment) the score for the tuple (word, t0, t1).
          auto &label_segments = accumulator_.emplace(
              label, MMT::mapped_type()).first->second;
          auto ret = label_segments.emplace(
              std::make_tuple(t0, t1), lkh_through_arc);
          if (!ret.second) {
            ret.first->second = LogAdd(ret.first->second, lkh_through_arc);
          }
        }
      }
    }

    elapsed_ = timer.Elapsed();
  }

 private:
  const std::string key_;
  CompactLattice clat_;
  const std::set<int32> exclude_labels_;
  BaseFloat acoustic_scale_, graph_scale_, insertion_penalty_, beam_;
  SegmentScoreWriter *score_writer_;
  std::vector<int32> state_times_;
  std::vector<double> fw_, bw_;
  MMT accumulator_;
  double total_lkh_;
  double elapsed_;
};
}  // namespace kaldi

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char *usage =
        "This tool is used to create a positional inverted index of the given "
        "lattices, where the score of each word in a segment is the "
        "probability that the word occurs in any of the transcriptions of the "
        "utterance at that specific time position (segment).\n"
        "\n"
        "Usage: lattice-word-index-segment [options] <lattice-rspecifier> "
        "<index-wspecifier>\n"
        " e.g.: lattice-word-index-segment --acoustic-scale=0.1 ark:1.lats "
        "ark:1.index\n";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    std::string exclude_symbols_str = "";
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                               "and adding the insertion penalty).");
    po.Register("exclude-words", &exclude_symbols_str,
                "Space-separated list of integers representing the words to "
                "exclude from the index.");
    // Register TaskSequencer options.
    TaskSequencerConfig task_sequencer_config;
    task_sequencer_config.Register(&po);

    po.Read(argc, argv);
    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    // Parse set of word symbols to exclude from the index
    set<kaldi::int32> exclude_symbols;
    {
      vector<kaldi::int32> tmp;
      kaldi::SplitStringToIntegers(exclude_symbols_str, " ", true, &tmp);
      kaldi::CopyVectorToSet(tmp, &exclude_symbols);
    }

    const std::string lattice_rspecifier = po.GetArg(1);
    LatticeScorerTask::SegmentScoreWriter score_writer(po.GetArg(2));

    TaskSequencer<LatticeScorerTask> task_sequencer(task_sequencer_config);
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      // Note: we will directly modify the lattice in the lattice_reader to
      // avoid copies.
      CompactLattice *clat = &lattice_reader.Value();
      // Schedule task.
      task_sequencer.Run(new LatticeScorerTask(lattice_key,
                                               *clat,
                                               exclude_symbols,
                                               acoustic_scale,
                                               graph_scale,
                                               insertion_penalty,
                                               beam,
                                               &score_writer));
      lattice_reader.FreeCurrent();
    }
    task_sequencer.Wait();
    score_writer.Close();
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

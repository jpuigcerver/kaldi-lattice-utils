// kwsbin2/lattice-word-index-position.cc

// Copyright (c) 2017-2018 Joan Puigcerver <joapuipe@upv.es>

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

#include "base/kaldi-common.h"
#include "base/timer.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils2.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"


namespace kaldi {
void ProcessLattice(
    const std::string &key, CompactLattice *clat, BaseFloat acoustic_scale,
    BaseFloat graph_scale, BaseFloat insertion_penalty, BaseFloat beam,
    std::vector<int32> *state_input_length, std::vector<int32> *state_time,
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

  Timer timer;
  if (clat->Start() != fst::kNoStateId) {
    // If needed, sort the compact lattice in topological order
    TopSortCompactLatticeIfNeeded(clat);
    // Make sure that all sequences arriving to each state have the same
    // length (number of words).
    CompactLattice tmp = *clat;
    DisambiguateStateInputSequenceLength(tmp, clat, state_input_length, false);
    *total_lkh = ComputeLatticeAlphasAndBetas(*clat, false, fw, bw);
    CompactLatticeStateTimes(*clat, state_time);
  } else {
    state_input_length->clear();
    fw->clear();
    bw->clear();
    *total_lkh = 0;
  }
  KALDI_VLOG(1) << "Lattice " << key << ": Preprocessing done in "
                << timer.Elapsed() << " seconds.";
}

class LatticeScorerTask {
 public:
  typedef TableWriter<
      BasicTupleVectorHolder<int32, int32, int32, int32, double>>
      PositionScoreWriter;
  typedef std::tuple<double, double, int32, int32> map_entry;
  typedef std::tuple<int32, int32, int32, int32, double> output_tuple;

  LatticeScorerTask(
      const std::string &key, const CompactLattice &clat,
      const std::set<int32> &include_labels,
      const std::set<int32> &exclude_labels,
      BaseFloat acoustic_scale, BaseFloat graph_scale,
      BaseFloat insertion_penalty, BaseFloat beam,
      PositionScoreWriter* posterior_writer) :
      key_(key), clat_(clat),
      include_labels_(include_labels),
      exclude_labels_(exclude_labels),
      acoustic_scale_(acoustic_scale), graph_scale_(graph_scale),
      insertion_penalty_(insertion_penalty), beam_(beam),
      score_writer_(posterior_writer) {}

  ~LatticeScorerTask() {
    // Put all positions of all words into a single vector.
    std::vector<output_tuple> position_scores;
    for (const auto& ws : accumulator_) {
      const int32 word_label = ws.first;
      const auto& positions = ws.second;
      for (const auto& pp : positions) {
        const auto pos = pp.first + 1;
        const double lkh = std::get<0>(pp.second);
        const auto t0 = std::get<2>(pp.second);
        const auto t1 = std::get<3>(pp.second);
        position_scores.emplace_back(word_label, pos, t0, t1, lkh - total_lkh_);
      }
    }
    // Sort the vector according to:
    //   1) Descending probability
    //   2) Ascending word label
    //   3) Ascending position
    std::sort(position_scores.begin(), position_scores.end(),
              [](const output_tuple &a, const output_tuple &b) -> bool {
                if (std::get<4>(b) != std::get<4>(a)) {
                  return std::get<4>(b) < std::get<4>(a);
                } else if (std::get<0>(a) != std::get<0>(b)) {
                  return std::get<0>(a) < std::get<0>(b);
                } else {
                  return std::get<1>(a) < std::get<1>(b);
                }
              });

    score_writer_->Write(key_, position_scores);

    KALDI_VLOG(1) << "Lattice " << key_ << ": done in " << elapsed_
                  << " seconds.";
  }

  void operator()() {
    // Prepare lattice to compute the scores.
    ProcessLattice(
        key_, &clat_, acoustic_scale_, graph_scale_, insertion_penalty_, beam_,
        &state_input_len_, &state_time_, &fw_, &bw_, &total_lkh_);

    // Timer does not include preprocessing time.
    Timer timer;

    for (fst::StateIterator<CompactLattice> sit(clat_); !sit.Done();
         sit.Next()) {
      for (fst::ArcIterator<CompactLattice> ait(clat_, sit.Value());
           !ait.Done(); ait.Next()) {
        const auto &arc = ait.Value();
        const auto &label = arc.olabel;
        const bool valid_label = (
            label != 0 && (
                (!include_labels_.empty() && include_labels_.count(label) > 0) ||
                    exclude_labels_.count(label) == 0
            )
        );
        if (valid_label) {
          // Word position within the sequence (0-based)
          const int32 pos = state_input_len_[sit.Value()];
          // Likelihood of the arc.
          const double arc_lkh = -ConvertToCost(arc.weight);
          // Likelihood of all paths through the arc.
          const double lkh_through_arc =
              fw_[sit.Value()] + arc_lkh + bw_[arc.nextstate];
          // Add (or increment) the score for the tuple (word, t0, t1).
          auto &label_positions = accumulator_.emplace(
              label, std::map<int32, map_entry>()).first->second;
          auto ret = label_positions.emplace(
              pos,
              std::make_tuple(lkh_through_arc,
                              lkh_through_arc,
                              state_time_[sit.Value()],
                              state_time_[arc.nextstate]));
          if (!ret.second) {
            const auto p = LogAdd(std::get<0>(ret.first->second), lkh_through_arc);
            auto a = std::get<1>(ret.first->second);
            auto t0 = std::get<2>(ret.first->second);
            auto t1 = std::get<3>(ret.first->second);
            if (lkh_through_arc > a) {
              a = lkh_through_arc;
              t0 = state_time_[sit.Value()];
              t1 = state_time_[arc.nextstate];
            }
            ret.first->second = std::make_tuple(p, a, t0, t1);
          }
        }
      }
    }

    elapsed_ = timer.Elapsed();
  }

 private:
  const std::string key_;
  CompactLattice clat_;
  const std::set<int32> include_labels_;
  const std::set<int32> exclude_labels_;
  BaseFloat acoustic_scale_, graph_scale_, insertion_penalty_, beam_;
  PositionScoreWriter *score_writer_;
  std::vector<int32> state_input_len_;
  std::vector<int32> state_time_;
  std::vector<double> fw_, bw_;
  std::map<int32, std::map<int32, map_entry>> accumulator_;
  double total_lkh_;
  double elapsed_;
};
}  // namespace kaldi

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage =
        "This tool creates a positional inverted index of the given lattices, "
        "in the traditional meaning of \"position\" in the context of search "
        "engines. That is, the probability that a word appears at some "
        "position within the transcription, for all possible transcriptions of "
        "the utterance.\n"
        "\n"
        "Usage: lattice-to-word-position-post [options] lat-rspecifier "
        "post-wspecifier\n"
        " e.g.: lattice-word-index-position --acoustic-scale=0.1 ark:1.lats "
        "ark:1.word.pos.post\n";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    std::string exclude_symbols_str, include_symbols_str;
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
    po.Register("include-words", &include_symbols_str,
                "Space-separated list of integers representing the words to "
                "include in the index.");
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

    // Parse set of word symbols to include in the index
    set<kaldi::int32> include_symbols;
    {
      vector<kaldi::int32> tmp;
      kaldi::SplitStringToIntegers(include_symbols_str, " ", true, &tmp);
      kaldi::CopyVectorToSet(tmp, &include_symbols);
    }

    const std::string lattice_rspecifier = po.GetArg(1);
    LatticeScorerTask::PositionScoreWriter score_writer(po.GetArg(2));

    TaskSequencer<LatticeScorerTask> task_sequencer(task_sequencer_config);
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      // Schedule task.
      task_sequencer.Run(new LatticeScorerTask(lattice_reader.Key(),
                                               lattice_reader.Value(),
                                               include_symbols,
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
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

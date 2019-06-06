// kwsbin2/lattice-char-index-segment.cc

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

#include "base/kaldi-common.h"
#include "base/timer.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils2.h"
#include "fstext/fst-info.h"
#include "kwsbin2/utils.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"

namespace kaldi {

void ProcessLattice(
    const std::string &key, CompactLattice *clat, BaseFloat acoustic_scale,
    BaseFloat graph_scale, BaseFloat insertion_penalty, BaseFloat beam) {
  Timer timer;
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
  }
  KALDI_VLOG(1) << "Lattice " << key << ": Preprocessing done in "
                << timer.Elapsed() << " seconds.";
}

template <typename Label, typename Int>
std::string LabelSequenceToStringAndPosition(
    const std::vector<Label>& ilabels,
    const std::vector<std::tuple<Label, Int>>& label_to_char_pos,
    bool skip_epsilon,
    Int* position) {
  KALDI_ASSERT(position != nullptr);
  *position = 0;

  // Convert label representing <char,word position> to char only,
  // and extract the position of the pseudo-word.
  std::vector<Label> tmp;
  for (const auto& ilabel : ilabels) {
    const auto olabel = std::get<0>(label_to_char_pos[ilabel]);
    const auto pos = std::get<1>(label_to_char_pos[ilabel]);
    tmp.push_back(olabel);
    if (*position == 0 && pos != 0) { *position = pos; }
  }

  return LabelSequenceToString(tmp, skip_epsilon);
}

template <typename Arc, typename Label, typename Int>
void PrintFst(
    const fst::Fst<Arc>& ifst,
    const std::vector<std::tuple<Label, Int>>& label_to_char_count,
    const std::vector<std::tuple<Int, Int>>& label_to_segm) {
  using namespace fst;
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto arc = ait.Value();
      std::cout << s << "\t"
                << arc.nextstate << "\t"
                << "(" << std::get<0>(label_to_char_count[arc.ilabel]) << ", "
                << std::get<1>(label_to_char_count[arc.ilabel]) << ")" << "\t"
                << "(" << std::get<0>(label_to_segm[arc.olabel]) << ", "
                << std::get<1>(label_to_segm[arc.olabel]) << ")" << "\t"
                << arc.weight.Value() << std::endl;
    }
  }
}

class LatticeScorerTask {
 public:
  typedef std::unordered_map<fst::StdArc::Label, fst::StdArc::Label> LabelMap;
  typedef std::unordered_set<fst::StdArc::Label> GroupSet;
  typedef std::tuple<std::string, int32, int32, int32, double> output_tuple;
  typedef TableWriter <
  BasicTupleVectorHolder<std::string, int32, int32, int32, double>>
      ScoreWriter;

  LatticeScorerTask(const std::string &key, const CompactLattice &clat,
                    const LabelMap& label_group,
                    const GroupSet& groups_inc_word_count,
                    const GroupSet& delete_groups,
                    const int32 num_nbest,
                    BaseFloat determinize_delta,
                    BaseFloat acoustic_scale, BaseFloat graph_scale,
                    BaseFloat insertion_penalty, BaseFloat beam,
                    ScoreWriter* score_writer) :
      key_(key), clat_(clat), label_group_(label_group),
      groups_inc_word_count_(groups_inc_word_count),
      delete_groups_(delete_groups), num_nbest_(num_nbest),
      determinize_delta_(determinize_delta), acoustic_scale_(acoustic_scale),
      graph_scale_(graph_scale), insertion_penalty_(insertion_penalty),
      beam_(beam), score_writer_(score_writer) {}

  ~LatticeScorerTask() {
    score_writer_->Write(key_, nbest_segment_);
    KALDI_VLOG(1) << "Lattice " << key_ << ": done in " << elapsed_
                  << " seconds.";
  }

  void operator()() {
    ProcessLattice(key_, &clat_, acoustic_scale_, graph_scale_,
                   insertion_penalty_, beam_);

    if (clat_.Start() == fst::kNoStateId) {
      elapsed_ = 0;
      return;
    }

    // Timer does not include preprocessing time.
    Timer timer;

    // Convert CompactLattice to LogFst with segmentation (t_start, t_end)
    // as the output label, and the (character, word position) as the input
    // label. Here, "word position" refers to the position within the sentence
    // of the word which the character belongs to.
    fst::VectorFst<fst::LogArc> tmp_fst;
    std::vector<fst::LogArc::Label> state_group;
    std::vector< std::tuple<fst::LogArc::Label, int32> > label_to_word_pos;
    std::vector< std::tuple<int32, int32> > label_to_segm;
    fst::CompactLatticeToWordCountSegmentFst(clat_,
                                             label_group_,
                                             groups_inc_word_count_,
                                             &tmp_fst,
                                             &state_group,
                                             &label_to_word_pos,
                                             &label_to_segm);
    clat_.DeleteStates();  // Free memory

    // Compute forward/backward
    fst::vector<fst::LogArc::Weight> fw, bw;
    fst::ShortestDistance(tmp_fst, &fw, /* reverse= */ false);
    fst::ShortestDistance(tmp_fst, &bw, /* reverse= */ true);
    const double total_cost = bw[tmp_fst.Start()].Value();

    // Create the FactorFst from the grouped states: A full-path in the
    // factor fst is equivalent to some sub-path in the original fst that
    // is made of only labels of the same group.
    fst::GroupFactorFst(&tmp_fst, state_group, fw, bw);

    // Remove arcs from non-interesting groups (e.g. whitespace).
    fst::DeleteArcs(
        &tmp_fst,
        // Lambda returns true if the label's group is in the delete_groups_.
        [&label_to_word_pos, this](const fst::LogArc& arc) -> bool {
          const auto label = std::get<0>(label_to_word_pos[arc.ilabel]);
          const auto it = label_group_.find(label);
          if (it != label_group_.end() && delete_groups_.count(it->second)) {
            return true;
          } else {
            return false;
          }
        });

    // We want word-segmentation, instead of the character-level segmentation.
    // Thus, we need to sum all the hypotheses with the same word-level
    // segmentation, even if the character-level segmentation is different.
    // In order to do so, we take advantage of the fact that every path from
    // The initial state to the final state is a word, thus we will remove
    // (set to epsilon) the output label of all arcs except those outgoing from
    // the initial state or entering a final state.
    fst::SymbolToPathSegmentationFst(&tmp_fst, label_to_segm);

    // Determinize in the log semiring to sum weights of all paths of the
    // same (word, position).
    fst::VectorFst<fst::LogArc> det_log;
    fst::Determinize(
        fst::ProjectFst<fst::LogArc>(tmp_fst, fst::PROJECT_INPUT),
        &det_log,
        fst::DeterminizeOptions<fst::LogArc>(
            determinize_delta_,
            fst::LogArc::Weight::Zero(),
            fst::kNoStateId,
            0,
            fst::DETERMINIZE_FUNCTIONAL));
    KALDI_VLOG(1) << "Lattice " << key_ << ": "
                  << fst::ComputeNumberOfPaths(&det_log)
                  << " pseudo-words.";

    // Determinize in the tropical semiring to obtain the best segmentation
    // for each (word, position).
    using AM1 = fst::WeightConvertMapper<fst::LogArc, fst::StdArc>;
    fst::VectorFst<fst::StdArc> det_std;
    fst::Determinize(
        fst::ArcMapFst<fst::LogArc, fst::StdArc, AM1>(tmp_fst, AM1()),
        &det_std,
        fst::DeterminizeOptions<fst::StdArc>(
            determinize_delta_,
            fst::StdArc::Weight::Zero(),
            fst::kNoStateId,
            0,
            fst::DETERMINIZE_DISAMBIGUATE));

    tmp_fst.DeleteStates();   // Free space.

    // Remove weights from the tropical fst
    fst::ArcMap(&det_std, fst::RmWeightMapper<fst::StdArc, fst::StdArc>());

    // Compose the two deterministic fsts to get a new deterministic fst with
    // the best segmentation and the total cost for each pseudo-word.
    fst::VectorFst<fst::StdArc> det_subpaths;
    fst::Compose(fst::ArcMapFst<fst::LogArc, fst::StdArc, AM1>(det_log, AM1()),
                 det_std, &det_subpaths);
    det_log.DeleteStates();   // Free space.
    det_std.DeleteStates();   // Free space.

    // Obtain n-best word segmentations
    fst::VectorFst<fst::StdArc> nbest_fst;
    fst::ShortestPath(det_subpaths, &nbest_fst, num_nbest_);

    // Process n-best list to prepare the output
    std::vector<fst::VectorFst<fst::StdArc>> nbest_fsts;
    ConvertNbestToVector(nbest_fst, &nbest_fsts);
    for (const auto& nbfst : nbest_fsts) {
      std::vector<fst::StdArc::Label> nb_isyms, nb_osyms;
      fst::StdArc::Weight nb_cost;
      GetLinearSymbolSequence(nbfst, &nb_isyms, &nb_osyms, &nb_cost);

      int32 position = 0;
      const auto pseudoword = LabelSequenceToStringAndPosition(
          nb_isyms, label_to_word_pos, /* skip_epsilon= */ true, &position);
      if (pseudoword.empty() || nb_osyms.empty()) {
        KALDI_WARN << "Lattice " << key_ << ": Ignoring eps pseudo-word";
        continue;
      }
      const auto t0 = nb_osyms.front() - 1;
      const auto t1 = nb_osyms.back() - 1;
      const auto prob = total_cost - nb_cost.Value();
      nbest_segment_.emplace_back(pseudoword, position, t0, t1, prob);
    }

    // Sort the vector according to:
    //   1) Descending probability
    //   2) Ascending word label
    //   3) Ascending word position
    std::sort(nbest_segment_.begin(), nbest_segment_.end(),
              [](const output_tuple &a, const output_tuple &b) -> bool {
                if (std::get<4>(a) != std::get<4>(b)) {
                  return std::get<4>(a) > std::get<4>(b);
                } else if (std::get<0>(a) != std::get<0>(b)) {
                  return std::get<0>(a) < std::get<0>(b);
                } else {
                  return std::get<1>(a) < std::get<1>(b);
                }
              });

    elapsed_ = timer.Elapsed();
  }

 private:
  const std::string key_;
  CompactLattice clat_;
  LabelMap label_group_;
  GroupSet groups_inc_word_count_;
  GroupSet delete_groups_;
  const int32 num_nbest_;
  BaseFloat determinize_delta_;
  BaseFloat acoustic_scale_, graph_scale_, insertion_penalty_, beam_;
  std::vector<output_tuple> nbest_segment_;
  ScoreWriter* score_writer_;
  double elapsed_;
};

}  // namespace kaldi


int main(int argc, char** argv) {
  try {
    using namespace kaldi;
    typedef fst::StdArc::Label Label;

    const char* usage =
        "Build a position-level word index from character lattices.\n"
        "\n"
        "Characters (arcs in the lattice) are grouped into groups, and "
        "subpaths from the lattice are extracted so that all arcs in the "
        "subpath are part of the same group.\n"
        "This creates a straightforward approach to obtain \"words\" from the "
        "character lattices, by grouping the whitespace characters and the "
        "rest of characters in two different groups.\n"
        "\n"
        "Characters (labels) that should not be indexed should be part of the "
        "\"wspace-group\" (e.g. different types of whitespaces).\n"
        "Other groups that should be considered as isolated words can also be "
        "specified using the \"--other-groups\" option. Labels within a group "
        "must be separated by spaces and groups are separated by a semicolon "
        "(e.g. --other-groups=\"1 2 3 ; 4 5 6 7\" creates two groups with 3 "
        "and 4 elements respectively). Labels not present in any specific "
        "group, are grouped together in a default group.\n"
        "\n"
        "Usage: lattice-char-index-segment [options] <wspace-group> "
        "<lattice-rspecifier> <index-wspecifier>\n"
        " e.g.: lattice-char-index-segment \"1 2\" ark:lats.ark ark,t:-\n"
        " e.g.: lattice-char-index-segment --nbest=10000 \"1 2\" "
        "ark:lats.ark ark:index.ark\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat determinize_delta = fst::kDelta / 8.0F;
    int nbest = 100;
    std::string other_groups_str = "";
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam,
                "Pruning beam (applied after acoustic scaling and adding the "
                "insertion penalty).");
    po.Register("nbest", &nbest, "Extract this number of n-best hypotheses.");
    po.Register("determinize-delta", &determinize_delta,
                "Delta threshold used for the determinization algorithm.");
    po.Register("other-groups", &other_groups_str,
                "Specific labels to group as words. Groups are separated "
                "with a semicolon, and labels within a group are separated "
                "by spaces.");
    // Register TaskSequencer options.
    TaskSequencerConfig task_sequencer_config;
    task_sequencer_config.Register(&po);

    po.Read(argc, argv);
    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    // Parse separator symbols from arguments
    std::unordered_set<Label> wspace_labels;
    std::unordered_map<Label, Label> label_group;
    std::unordered_set<Label> groups_inc_word_count;
    ParseSeparatorGroups(po.GetArg(1), other_groups_str, &wspace_labels,
                         &label_group, &groups_inc_word_count);

    Timer timer;
    int num_lats = 0;
    LatticeScorerTask::ScoreWriter score_writer(po.GetArg(3));
    TaskSequencer<LatticeScorerTask> task_sequencer(task_sequencer_config);
    for (SequentialCompactLatticeReader lattice_reader(po.GetArg(2));
         !lattice_reader.Done(); lattice_reader.Next(), ++num_lats) {
      const std::string lattice_key = lattice_reader.Key();
      // Schedule task.
      task_sequencer.Run(new LatticeScorerTask(lattice_reader.Key(),
                                               lattice_reader.Value(),
                                               label_group,
                                               groups_inc_word_count,
                                               std::unordered_set<Label>{1},
                                               nbest,
                                               determinize_delta,
                                               acoustic_scale,
                                               graph_scale,
                                               insertion_penalty,
                                               beam,
                                               &score_writer));
      lattice_reader.FreeCurrent();
    }
    task_sequencer.Wait();
    score_writer.Close();

    const auto effective_num_threads =
        std::min<size_t>(num_lats, task_sequencer_config.num_threads);
    KALDI_LOG << "Done " << num_lats << " lattices in "  << timer.Elapsed()
              << "s using " << effective_num_threads << " threads.";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}

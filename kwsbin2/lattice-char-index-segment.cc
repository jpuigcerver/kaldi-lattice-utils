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

template <typename Container>
std::string LabelSequenceToString(const Container& label_sequence,
                                  bool skip_epsilon) {
  size_t i = 0;
  // Skip <eps> symbols at the start of the word.
  if (skip_epsilon) {
    for (; i < label_sequence.size() && label_sequence[i] == 0; ++i) {}
    if (i == label_sequence.size()) return "";
  }
  std::string pseudoword = std::to_string(label_sequence[i]);
  for (size_t i = 1; i < label_sequence.size(); ++i) {
    if (skip_epsilon && label_sequence[i] == 0) continue;
    pseudoword += "_" + std::to_string(label_sequence[i]);
  }
  return pseudoword;
}

void ProcessLattice(
    const std::string &key, CompactLattice *clat, BaseFloat acoustic_scale,
    BaseFloat graph_scale, BaseFloat insertion_penalty, BaseFloat beam) {
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
}

class LatticeScorerTask {
 public:
  typedef std::unordered_map<fst::StdArc::Label, fst::StdArc::Label> LabelMap;
  typedef std::unordered_set<fst::StdArc::Label> GroupSet;
  typedef TableWriter<BasicTupleVectorHolder<std::string, double, int32, int32>>
      SegmentScoreWriter;

  LatticeScorerTask(const std::string &key, const CompactLattice &clat,
                    BaseFloat acoustic_scale, BaseFloat graph_scale,
                    BaseFloat insertion_penalty, BaseFloat beam,
                    const LabelMap& label_group,
                    const GroupSet& delete_groups,
                    const int32 num_nbest,
                    SegmentScoreWriter* score_writer) :
      key_(key), clat_(clat), acoustic_scale_(acoustic_scale),
      graph_scale_(graph_scale), insertion_penalty_(insertion_penalty),
      beam_(beam), label_group_(label_group), delete_groups_(delete_groups),
      num_nbest_(num_nbest), score_writer_(score_writer) {}

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
    // as the output label, and the characters as the input label.
    fst::VectorFst<fst::LogArc> tmp_fst, subpath_fst;
    std::vector<std::tuple<int32, int32>> label_to_segm;
    CompactLatticeToSegmentFst(clat_, &tmp_fst, &label_to_segm);
    clat_.DeleteStates();  // Free memory

    // Disambiguate states in the fst, so that all input labels to a
    // state are part of the same group.
    std::vector<fst::LogArc::Label> tmp_state_group;
    fst::DisambiguateStatesByInputLabelGroup(
        tmp_fst, label_group_, &subpath_fst, &tmp_state_group,
        /* use_input= */ true);
    tmp_fst.DeleteStates();  // Free memory

    // Compute forward/backward
    fst::vector<fst::LogArc::Weight> fw, bw;
    fst::ShortestDistance(subpath_fst, &fw, /* reverse= */ false);
    fst::ShortestDistance(subpath_fst, &bw, /* reverse= */ true);
    const double total_cost = bw[subpath_fst.Start()].Value();

    // Create the FactorFst from the grouped states: A full-path in the
    // factor fst is equivalent to some sub-path in the original fst that
    // is made of only labels of the same group.
    fst::GroupFactorFst(&subpath_fst, tmp_state_group, fw, bw);

    // Remove arcs from non-interesting groups (e.g. whitespace).
    fst::RemoveArcsByGroup(
        &subpath_fst, label_group_, delete_groups_, /* use_input= */ true);

    // We want word-segmentation, instead of the character-level segmentation.
    // Thus, we need to sum all the hypotheses with the same word-level
    // segmentation, even if the character-level segmentation is different.
    // In order to do so, we take advantage of the fact that every path from
    // The initial state to the final state is a word, thus we will remove
    // (set to epsilon) the output label of all arcs except those outgoing from
    // the initial state or entering a final state.
    fst::SymbolToPathSegmentationFst(&subpath_fst, label_to_segm);

    // We need to sum up all scores for the same word segmentation.
    // That means determinization in the log-semiring. However, since
    // the FST is non-functional we need to encode the (ilabel, olabel)
    // pairs into a new label and make the original FST a weighted
    // automaton.
    fst::EncodeMapper<fst::LogArc> encoder(fst::kEncodeLabels, fst::ENCODE);
    fst::Encode(&subpath_fst, &encoder);

    // Determinize subpaths FST
    fst::VectorFst<fst::LogArc> det_subpath_fst;
    fst::Determinize(subpath_fst, &det_subpath_fst);
    subpath_fst.DeleteStates();  // Free memory

    KALDI_VLOG(1) << "Lattice " << key_ << ": "
                  << fst::ComputeNumberOfPaths(&det_subpath_fst)
                  << " pseudo-words.";

    // Recover segmentation information.
    fst::Decode(&det_subpath_fst, encoder);

    // Convert Fst from log-semiring to tropical-semiring
    fst::VectorFst<fst::StdArc> det_subpath_fst_std;
    fst::ArcMap(det_subpath_fst, &det_subpath_fst_std,
                fst::WeightConvertMapper<fst::LogArc, fst::StdArc>());
    det_subpath_fst.DeleteStates();  // Free memory

    // Obtain n-best word segmentations
    fst::VectorFst<fst::StdArc> nbest_fst;
    fst::ShortestPath(det_subpath_fst_std, &nbest_fst, num_nbest_);

    // Process n-best list to prepare the output
    std::vector<fst::VectorFst<fst::StdArc>> nbest_fsts;
    ConvertNbestToVector(nbest_fst, &nbest_fsts);
    for (const auto& nbfst : nbest_fsts) {
      std::vector<fst::StdArc::Label> nb_isyms, nb_osyms;
      fst::StdArc::Weight nb_cost;
      GetLinearSymbolSequence(nbfst, &nb_isyms, &nb_osyms, &nb_cost);

      const auto pseudoword =
          LabelSequenceToString(nb_isyms, /* skip_epsilon= */ true);
      if (pseudoword.empty() || nb_osyms.empty()) {
        KALDI_WARN << "Lattice " << key_ << ": Ignoring eps pseudo-word";
        continue;
      }
      const auto t0 = nb_osyms.front() - 1;
      const auto t1 = nb_osyms.back() - 1;
      const auto prob = total_cost - nb_cost.Value();
      nbest_segment_.emplace_back(pseudoword, prob, t0, t1);
    }

    elapsed_ = timer.Elapsed();
  }

 private:
  const std::string key_;
  CompactLattice clat_;
  BaseFloat acoustic_scale_, graph_scale_, insertion_penalty_, beam_;
  LabelMap label_group_;
  GroupSet delete_groups_;
  std::vector<std::tuple<std::string, double, int32, int32>> nbest_segment_;
  const int32 num_nbest_;
  SegmentScoreWriter* score_writer_;
  double elapsed_;
};

}  // namespace kaldi


int main(int argc, char** argv) {
  try {
    using namespace kaldi;
    typedef fst::StdArc::Label Label;

    const char* usage =
        "Build a word index from character lattices.\n"
            "\n"
            "Words are any sequence of characters in between any of the separator "
            "symbols. The program will output the n-best character segmentations "
            "of words, with their score. More precisely:\n"
            "Let R_c = 1 denote y \\in (\\Sigma^* @)? c_{1:n} (@ \\Sigma^*)? and "
            "R_c = 0 otherwise, where @ is any of the separator symbols.\n"
            "Then R denotes whether the transcription (y) of each utterance "
            "contains the word formed by characters c_1,...,c_n.\n"
            "Let s_{1:n} be a segmentation of each character in c_{1:n}, then "
            "the program computes:\n"
            "\n"
            "1. If --only-best-segmentation=false (the default) then:\n"
            "n-best_{c_{1:n},s_{1:n}} P(R_c = 1| x, s_{1:n})\n"
            "\n"
            "2. If --only-best-segmentation=true then:\n"
            "n-best_{c_{1:n},s_{1:n}} max_{s_{1:n}} P(R_c = 1 | x, s_{1:n})\n"
            "\n"
            "This gives a lower bound to P(R_c = 1 | x), but it is usually quite "
            "close.\n"
            "\n"
            "Usage: kaldi-lattice-word-index [options] separator-symbols "
            "lat-rspecifier\n"
            " e.g.: kaldi-lattice-word-index \"1 2\" ark:lats.ark\n"
            " e.g.: kaldi-lattice-word-index --nbest=10000 \"1 2\" ark:lats.ark\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    int nbest = 100;

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
    ParseSeparatorGroups(po.GetArg(1), po.GetArg(2), &wspace_labels,
                         &label_group, &groups_inc_word_count);

    LatticeScorerTask::SegmentScoreWriter score_writer(po.GetArg(2));
    TaskSequencer<LatticeScorerTask> task_sequencer(task_sequencer_config);
    for (SequentialCompactLatticeReader lattice_reader(po.GetArg(3));
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      // Schedule task.
      task_sequencer.Run(new LatticeScorerTask(lattice_reader.Key(),
                                               lattice_reader.Value(),
                                               acoustic_scale,
                                               graph_scale,
                                               insertion_penalty,
                                               beam,
                                               label_group,
                                               std::unordered_set<Label>{1},
                                               nbest,
                                               &score_writer));
      lattice_reader.FreeCurrent();
    }
    task_sequencer.Wait();
    score_writer.Close();
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}

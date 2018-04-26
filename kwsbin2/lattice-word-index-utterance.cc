// kwsbin2/lattice-word-index-utterance.cc

// Copyright (c) 2017 Joan Puigcerver <joapuipe@upv.es>

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
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"
#include "util/common-utils.h"
#include "util/kaldi-thread.h"

namespace fst {

template<typename Arc>
void CreateQueryFst(const typename Arc::Label &query_label,
                    const typename Arc::Label &rho_label,
                    MutableFst<Arc> *fst) {
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;

  fst->DeleteStates();
  fst->AddState();
  fst->AddState();
  fst->SetStart(0);
  fst->SetFinal(0, Weight::Zero());
  fst->SetFinal(1, Weight::One());

  // Add arcs sorted.
  if (query_label < rho_label) {
    fst->AddArc(0, Arc(query_label, query_label, Weight::One(), 1));
    fst->AddArc(0, Arc(rho_label, rho_label, Weight::One(), 0));
  } else {
    fst->AddArc(0, Arc(rho_label, rho_label, Weight::One(), 0));
    fst->AddArc(0, Arc(query_label, query_label, Weight::One(), 1));
  }
  fst->AddArc(1, Arc(rho_label, rho_label, Weight::One(), 1));

  fst->SetProperties(kFstProperties, fst->Properties(kFstProperties, true));
}

}  // namespace fst

namespace kaldi {

template<typename Set>
void FilterSymbols(const std::vector<int32> &symbols,
                   const Set &exclude_symbols,
                   std::vector<int32> *filtered_symbols) {
  filtered_symbols->reserve(symbols.size());
  for (const auto &s : symbols) {
    if (exclude_symbols.find(s) == exclude_symbols.end()) {
      filtered_symbols->push_back(s);
    }
  }
}

void ProcessLattice(
    const std::string &key, CompactLattice *clat, BaseFloat acoustic_scale,
    BaseFloat graph_scale, BaseFloat insertion_penalty, BaseFloat beam,
    const std::set<int32> &excluded_word_symbols, double *total_lkh,
    std::vector<int32> *all_word_symbols,
    std::vector<int32> *interesting_word_symbols) {
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

  all_word_symbols->clear();
  interesting_word_symbols->clear();
  if (clat->Start() != fst::kNoStateId) {
    // If needed, sort the compact lattice in topological order
    TopSortCompactLatticeIfNeeded(clat);
    // Compute lattice total likelihood
    std::vector<double> bw_lkh;
    ComputeCompactLatticeBetas(*clat, &bw_lkh);
    *total_lkh = bw_lkh[clat->Start()];
    GetOutputSymbols(*clat, false /* do not include epsilon */,
                     all_word_symbols);
    FilterSymbols(*all_word_symbols, excluded_word_symbols,
                  interesting_word_symbols);
  } else {
    *total_lkh = kLogZeroDouble;
  }
}

class ComputeWordScoreTask {
 public:

  ComputeWordScoreTask(
      const std::string &lattice_key,
      const CompactLattice &clat,
      const int32 word_symbol,
      const int32 rho_symbol,
      const double total_lkh,
      std::vector<std::tuple<int32, double>> *utterance_scores)
      : key_(lattice_key), clat_(clat), word_symbol_(word_symbol),
        rho_symbol_(rho_symbol), total_lkh_(total_lkh),
        utterance_scores_(utterance_scores) {}

  ~ComputeWordScoreTask() {
    if (query_lkh_ > total_lkh_) {
      KALDI_WARN << "Lattice " << key_ << ": found word (" << word_symbol_
                 << ") with higher likelihood than the total likelihood ("
                 << query_lkh_ << " vs. " << total_lkh_ << ")";
    }
    utterance_scores_->emplace_back(word_symbol_, query_lkh_ - total_lkh_);
  }

  void operator()() {
    CompactLattice query_lat, composed_lat;
    // Create the fst representing the language \Sigma^* symbol \Sigma^*
    // The query fst must be deterministic so that no repeated paths are
    // added in the composition with the log_fst (i.e. original lattice).
    CreateQueryFst(word_symbol_, rho_symbol_, &query_lat);
    // Compose / intersect with the lattice with the query fst to obtain
    // all paths in the lattice containing the query word.
    fst::RhoCompose(clat_, query_lat, rho_symbol_, &composed_lat);
    if (composed_lat.Start() != fst::kNoStateId) {
      // Sort states in topological order, if needed.
      TopSortCompactLatticeIfNeeded(&composed_lat);
      // Compute total likelihood of the intersection.
      std::vector<double> bw_lkh;
      ComputeCompactLatticeBetas(composed_lat, &bw_lkh);
      query_lkh_ = bw_lkh[composed_lat.Start()];
    } else {
      query_lkh_ = kLogZeroDouble;
    }
  }

 private:
  const std::string key_;
  const CompactLattice &clat_;
  int32 word_symbol_;
  int32 rho_symbol_;
  double total_lkh_;
  double query_lkh_;
  std::vector<std::tuple<int32, double>> *utterance_scores_;
  double elapsed_;
};

}  // namespace kaldi

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char *usage =
        "Create an inverted index of the given lattices, where the score of "
        "each word is the probability that the word is in the transcription "
        "of the signal. I.e. P(R = 1 | x, v), where x is the image and v is "
        "the word.\n"
        "\n"
        "For additional details, check the paper \"Probabilistic "
        "interpretation and improvements to the HMM-filler for handwritten "
        "keyword spotting\", by J. Puigcerver et al.\n"
        "\n"
        "Usage: lattice-word-index-utterance [options] <lattice-rspecifier> "
        "<index-wspecifier>\n"
        " e.g.: lattice-word-index-utterance --acoustic-scale=0.1 ark:1.lats "
        "ark:1.index\n";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    std::string exclude_symbols_str = "";
    int32 rho_symbol = std::numeric_limits<int32>::max();
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
    po.Register("rho-label", &rho_symbol,
                "Label that represents all possible word labels, it must not "
                "be equal to any other label (you shouldn't need to modify this "
                "in normal cases).");

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

    typedef TableWriter<BasicTupleVectorHolder<int32, double>>
        UtteranceScoreWriter;

    const std::string lattice_rspecifier = po.GetArg(1);
    UtteranceScoreWriter utterance_writer(po.GetArg(2));

    TaskSequencer<ComputeWordScoreTask> task_sequencer(task_sequencer_config);
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      // Note: we will directly modify the lattice in the lattice_reader to
      // avoid copies.
      CompactLattice *clat = &lattice_reader.Value();
      // Prepare lattice to compute the score of all queries.
      double total_lkh;
      std::vector<int32> all_word_symbols, interesting_word_symbols;
      ProcessLattice(
          lattice_key, clat, acoustic_scale, graph_scale, insertion_penalty,
          beam, exclude_symbols, &total_lkh, &all_word_symbols,
          &interesting_word_symbols);
      // Schedule tasks to compute the scores of each word symbol.
      std::vector<std::tuple<int32, double>> utterance_scores;
      Timer timer;
      for (const auto &word_symbol : interesting_word_symbols) {
        task_sequencer.Run(new ComputeWordScoreTask(
            lattice_key,
            *clat,
            word_symbol,
            rho_symbol,
            total_lkh,
            &utterance_scores));
      }
      task_sequencer.Wait();
      // Sort scores in decreasing order.
      std::sort(utterance_scores.begin(), utterance_scores.end(),
                [](const std::tuple<int32, double> &a,
                   const std::tuple<int32, double> &b) -> bool {
                  if (std::get<1>(b) != std::get<1>(a)) {
                    return std::get<1>(b) < std::get<1>(a);
                  } else {
                    return std::get<0>(a) < std::get<0>(b);
                  }
                });
      // Write scores to the table
      utterance_writer.Write(lattice_key, utterance_scores);
      lattice_reader.FreeCurrent();
      // Print running time stats.
      const auto effective_num_threads =
          std::min<size_t>(interesting_word_symbols.size(),
                   task_sequencer_config.num_threads);
      KALDI_VLOG(1) << "Lattice " << lattice_key << ": done in "
                    << timer.Elapsed() << " seconds for "
                    << interesting_word_symbols.size() << " words "
                    << "using " << effective_num_threads << " num threads.";
    }
    utterance_writer.Close();
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
}

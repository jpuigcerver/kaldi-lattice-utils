// kwsbin2/lattice-index.cc

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
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"
#include "util/basic-tuple-vector-holder.h"

#include "fst/script/print.h"

namespace fst {

void ConvertLatticeWeight(const kaldi::CompactLatticeWeight& iw,
                          Log64Weight* ow) {
  KALDI_ASSERT(ow != NULL);
  const Log64Weight w1(iw.Weight().Value1());
  const Log64Weight w2(iw.Weight().Value2());
  *ow = Times(w1, w2);
}

template <typename Arc>
void GetSymbols(const Fst<Arc>& fst,
                std::set<typename Arc::Label>* symbols,
                bool input = false) {
  symbols->clear();
  for (StateIterator< Fst<Arc> > siter(fst); !siter.Done(); siter.Next()) {
    for (ArcIterator< Fst<Arc> > aiter(fst, siter.Value()); !aiter.Done();
         aiter.Next()) {
      const Arc& arc = aiter.Value();
      const typename Arc::Label label = input ? arc.ilabel : arc.olabel;
      if (label != 0) { symbols->insert(label); }
    }
  }
}

template <typename Arc>
void CreateQueryFst(const std::set<typename Arc::Label>& alphabet,
                    const typename Arc::Label& query_label,
                    MutableFst<Arc>* fst) {
  assert(alphabet.count(query_label) != 0);
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;

  fst->DeleteStates();
  fst->AddState();
  fst->AddState();
  fst->SetStart(0);

  for (const Label& label : alphabet) {
    if (label == 0) continue;
    if (label == query_label) {
      fst->AddArc(0, Arc(label, label, Weight::One(), 1));
    } else {
      fst->AddArc(0, Arc(label, label, Weight::One(), 0));
    }
    fst->AddArc(1, Arc(label, label, Weight::One(), 1));
  }

  fst->SetFinal(0, Weight::Zero());
  fst->SetFinal(1, Weight::One());

  constexpr uint64_t properties =
      kExpanded | kMutable | kAcceptor | kIDeterministic | kODeterministic |
      kNoEpsilons | kNoIEpsilons | kNoOEpsilons | kUnweighted |
      kInitialCyclic | kAccessible | kCoAccessible;
  fst->SetProperties(kFstProperties, properties);
}

template <typename Arc>
int32 CompactLatticeToSegmentFst(const kaldi::CompactLattice& clat,
                                 MutableFst<Arc>* fst,
                                 std::vector<tuple<int32, int32>>* label2segm) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  using kaldi::CompactLattice;
  using kaldi::CompactLatticeWeight;
  using kaldi::LatticeWeight;
  // Assumption: The lattice must be aligned!
  // Compute the times of each state.
  std::vector<int32> times;
  const int32 total_frames = kaldi::CompactLatticeStateTimes(clat, &times);

  // TODO(jpuigcerver): Switch to unordered_map ?
  std::map<std::tuple<int32, int32>, Label> segm2label;

  for (StateIterator<CompactLattice> siter(clat); !siter.Done(); siter.Next()) {
    const StateId u = siter.Value();
    for (ArcIterator<CompactLattice> aiter(clat, u); !aiter.Done();
         aiter.Next()) {
      const auto& arc = aiter.Value();
      const StateId v = arc.nextstate;
      const Weight w = Weight(arc.weight.Weight().Value1() +
                              arc.weight.Weight().Value2());
      const auto segment = std::make_tuple(times[u], times[v]);
      const Label ilabel = segm2label.emplace(segment, segm2label.size() + 1).first->second;
      const Label olabel = arc.olabel;
      fst->AddArc(u, Log64Arc(ilabel, olabel, w, v));
    }

    if (clat.Final(u) != CompactLatticeWeight::Zero()) {
      const LatticeWeight& cw = clat.Final(u).Weight();
      const Weight w = Weight(cw.Value1() + cw.Value2());
      fst->SetFinal(u, w);
    }
  }

  // Reverse the map into a vector.
  // segm2label table starts from 1, thus the +1 here.
  label2segm->resize(segm2label.size() + 1);
  for (const auto& segm_label : segm2label) {
    (*label2segm)[segm_label.second] = segm_label.first;
  }

  return total_frames;
}

}  // namespace fst

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using namespace fst;

    const char* usage = "";

    ParseOptions po(usage);
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    double error_threshold = 0.01;
    std::string index_type = "utterance";
    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("index-type", &index_type,
                "Specify the type of the index to build: \"utterance\", "
                "\"segment\" or \"frame\".");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Read(argc, argv);

    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    // Index type: utterance-level, segment-level and frame-level
    std::transform(index_type.begin(), index_type.end(), index_type.begin(),
                   ::tolower);
    if (index_type != "utterance" && index_type != "segment" &&
        index_type != "frame") {
      KALDI_ERR << "Unknown --index-type \"" << index_type << "\"!";
      return -1;
    }

    typedef TableWriter<BasicTupleVectorHolder<int32, double>>
        UtteranceScoreWriter;
    typedef TableWriter<BasicTupleVectorHolder<double, int32, int32, int32>>
        SegmentScoreWriter;

    UtteranceScoreWriter* utterance_writer = nullptr;
    SegmentScoreWriter* segment_writer = nullptr;
    PosteriorWriter* posterior_writer = nullptr;
    if (index_type == "utterance") {
      utterance_writer = new UtteranceScoreWriter(po.GetArg(2));
    } else if (index_type == "segment") {
      segment_writer = new SegmentScoreWriter(po.GetArg(2));
    } else if (index_type == "frame") {
      posterior_writer = new PosteriorWriter(po.GetArg(2));
    }

    const auto lattice_scale = LatticeScale(graph_scale, acoustic_scale);

    const std::string lattice_rspecifier = po.GetArg(1);
    for (SequentialCompactLatticeReader lattice_reader(lattice_rspecifier);
         !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();
      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        ScaleLattice(lattice_scale, &clat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddWordInsPenToCompactLattice(insertion_penalty, &clat);
      // Lattice prunning
      if (beam != std::numeric_limits<BaseFloat>::infinity())
        PruneLattice(beam, &clat);
      // If needed, sort the compact lattice in topological order
      TopSortCompactLatticeIfNeeded(&clat);
      // DON'T DO THIS! IT CAN BREAK ALIGNMENT WHEN EPSILON SYMBOLS ADDED BY
      // lattice-align-words-lexicon:
      // If needed, remove epsilons from the compact lattice
      //if (clat.Properties(kEpsilons, true) == kEpsilons) {
      //  RmEpsilon(&clat);
      //}

      if (index_type == "utterance") {
        // Convert CompactLattice to Fst in the LogSemiring.
        // Note: This removes the alignment implicity.
        VectorFst<Log64Arc> log_fst;
        ConvertLattice(clat, &log_fst);
        // Compute total backward cost of each state, in order to compute the
        // total cost in the lattice (cost = - log-likelihood).
        std::vector<Log64Weight> distance;
        ShortestDistance(log_fst, &distance, true);
        const double total_cost = distance[log_fst.Start()].Value();
        // Get all the words in the lattice
        std::set<Log64Arc::Label> symbols;
        GetSymbols(log_fst, &symbols, false /* Get output symbols */);

        std::vector<std::tuple<int32, double>> utterance_scores;
        for (const Log64Arc::Label& symbol : symbols) {
          VectorFst<Log64Arc> qfst, cfst;
          // Create the fst representing the language \Sigma^* symbol \Sigma^*
          // The query fst must be deterministic so that no repeated paths are
          // added in the composition with the log_fst (i.e. original lattice).
          CreateQueryFst(symbols, symbol, &qfst);
          // Compose / intersect with the log_fst to obtain all paths in the
          // language represented by the query fst.
          Compose(log_fst, qfst, &cfst);
          // Compute the total cost in the resulting fst, which is the cost of the
          // paths accepting the query fst.
          ShortestDistance(cfst, &distance, true);
          const double query_cost = distance[cfst.Start()].Value();
          if (total_cost - query_cost > error_threshold) {
            KALDI_WARN << "Found query with higher likelihood than the "
                       << "likelihood (" << -query_cost << " vs. "
                       << -total_cost << ")";
          }
          utterance_scores.emplace_back(symbol, total_cost - query_cost);
        }

        // Sort scores in decreasing order.
        std::sort(utterance_scores.begin(), utterance_scores.end(),
                  [](const std::tuple<int32, double>& a,
                     const std::tuple<int32, double>& b) -> bool {
                    if (std::get<1>(b) != std::get<1>(a)) {
                      return std::get<1>(b) < std::get<1>(a);
                    } else {
                      return std::get<0>(a) < std::get<0>(b);
                    }
                  });
        // Write scores to the table
        utterance_writer->Write(lattice_key, utterance_scores);

      } else if (index_type == "segment") {
        // Convert the CompactLattice into a fst where the input labels
        // represent the segment in which the word appeared in the utterance.
        // Cost: O(V + E)
        VectorFst<Log64Arc> log_fst;
        std::vector<std::tuple<int32, int32>> label2segm;
        CompactLatticeToSegmentFst(clat, &log_fst, &label2segm);
        // Compute forward and backward of each state.
        // Cost: O(V + E)
        std::vector<Log64Weight> fw_cost, bw_cost;
        ShortestDistance(log_fst, &fw_cost, false);
        ShortestDistance(log_fst, &bw_cost, true);
        const double total_cost = bw_cost[log_fst.Start()].Value();

        map<std::tuple<Log64Arc::Label, Log64Arc::Label>, Log64Weight> word_segm_cost;
        for (StateIterator<VectorFst<Log64Arc>> siter(log_fst); !siter.Done();
             siter.Next()) {
          const auto u = siter.Value();
          for (ArcIterator<VectorFst<Log64Arc>> aiter(log_fst, u); !aiter.Done();
               aiter.Next()) {
            const auto& arc = aiter.Value();
            const auto ilabel = arc.ilabel;
            const auto olabel = arc.olabel;
            const auto v = arc.nextstate;
            // Total cost through the arc.
            const auto w = Times(arc.weight, Times(fw_cost[u], bw_cost[v]));
            auto r =
                word_segm_cost.emplace(std::make_tuple(ilabel, olabel), w);
            if (!r.second) {
              r.first->second = Plus(r.first->second, w);
            }
          }
        }

        typedef std::tuple<double, Log64Arc::Label, int32, int32> RTuple;
        std::vector<RTuple> cost_word_segm;
        cost_word_segm.reserve(word_segm_cost.size());
        for (const auto& kv : word_segm_cost) {
          const auto& ilabel = std::get<0>(kv.first);
          const auto& olabel = std::get<1>(kv.first);
          const int32 fbeg = std::get<0>(label2segm[ilabel]);
          const int32 fend = std::get<1>(label2segm[ilabel]);
          const double cost = kv.second.Value();
          cost_word_segm.push_back(std::make_tuple(cost, olabel, fbeg, fend));
        }
        std::sort(cost_word_segm.begin(), cost_word_segm.end(),
                  std::greater<RTuple>());

        segment_writer->Write(lattice_key, cost_word_segm);

      } else if (index_type == "frame") {
        // Compute the times of each state.
        // Assumption: The lattice must be aligned!
        std::vector<int32> times;
        const int32 total_frames = kaldi::CompactLatticeStateTimes(clat, &times);
        // Convert CompactLattice to a fst in the log semiring.
        // Note: This removes the alignment implicity.
        VectorFst<Log64Arc> log_fst;
        ConvertLattice(clat, &log_fst);

        bool misaligned_lattice = false;
        std::vector<std::unordered_map<Log64Arc::Label, Log64Weight>> accumulator(total_frames);
        for (StateIterator<VectorFst<Log64Arc>> siter(log_fst); !siter.Done();
             siter.Next()) {
          const auto u = siter.Value();
          for (ArcIterator<VectorFst<Log64Arc>> aiter(log_fst, u); !aiter.Done();
               aiter.Next()) {
            const auto& arc = aiter.Value();
            const auto olabel = arc.olabel;
            const auto v = arc.nextstate;
            const auto w = arc.weight;
            if (times[u] == times[v] && !misaligned_lattice) {
              misaligned_lattice = true;
              KALDI_WARN << "Lattice " << lattice_key << " is misaligned, "
                         << "a word with zero duration was found!";
            }
            for (int32 k = times[u]; k < times[v]; ++k) {
              auto r = accumulator[k].emplace(olabel, w);
              if (!r.second) {
                r.first->second = Plus(r.first->second, w);
              }
            }
          }
        }

        // Accumulate total cost through each frame
        std::vector<Log64Weight> total_cost_frame(total_frames, Log64Weight::Zero());
        for (int32 k = 0; k < total_frames; ++k) {
          for (const auto& kv : accumulator[k]) {
            total_cost_frame[k] = Plus(total_cost_frame[k], kv.second);
          }
        }

        PosteriorWriter::T posterior(total_frames);
        for (int32 k = 0; k < total_frames; ++k) {
          for (const auto& kv : accumulator[k]) {
            const int32 label = kv.first;
            const float logp  = total_cost_frame[k].Value() - kv.second.Value();
            posterior[k].emplace_back(label, logp);
          }
          std::sort(posterior[k].begin(), posterior[k].end());
        }

        posterior_writer->Write(lattice_key, posterior);
      }
    }

    if (utterance_writer != nullptr) {
      utterance_writer->Close();
      delete utterance_writer;
    } else if (segment_writer != nullptr) {
      segment_writer->Close();
      delete segment_writer;
    } else if (posterior_writer != nullptr) {
      posterior_writer->Close();
      delete posterior_writer;
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}

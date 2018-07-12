#undef NDEBUG

#include "fstext/expand-subpaths-between-delimiters.h"
#include "fstext/expand-subpaths-test-utils.h"
#include "fstext/fst-test-utils2.h"
#include "fstext/rand-fst.h"
#include "fstext/kaldi-fst-io.h"

using namespace fst;

template <typename Arc>
void RandFstSpecial(
    VectorFst<Arc>* ofst, std::set<typename Arc::Label>* delimiters) {
  typedef typename Arc::Weight Weight;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Label Label;
  ofst->DeleteStates();
  delimiters->clear();

  const int num_states = 3 + rand() % 8;  // values in [3, 4, ..., 10]
  for (int s = 0; s < num_states; ++s) {
    ofst->AddState();
    // 10% of states are final
    if (rand() % 10 == 0) {
      ofst->SetFinal(s, Weight((rand() % 5) / 4.0));
    }
  }
  // First state is always initial
  ofst->SetStart(0);
  // Last state is always final
  ofst->SetFinal(num_states - 1, Weight(0.0));

  const int num_total_labels = 3 + rand() % 5;  // values in [3, 4, 5, 6, 7]

  // Choose delimiters
  const int num_delimiters = rand() % num_total_labels;  // values [0...L-1]
  while (delimiters->size() < num_delimiters) {
    delimiters->emplace(1 + rand() % (num_total_labels - 1));
  }

  // Choose non-delimiters (the rest of labels, except epsilon)
  std::set<typename Arc::Label> regular;
  for (int label = 1; label < num_total_labels; ++label) {
    if (delimiters->count(label) == 0) {
      regular.emplace(label);
    }
  }

  std::set<StateId> states_with_input_delim;
  std::set<StateId> states_with_input_regular;
  for (StateId s = 1; s < num_states; ++s) {
    if(rand() % 2 == 0) states_with_input_delim.emplace(s);
    if(rand() % 2 == 0) states_with_input_regular.emplace(s);
  }

  auto add_arc_to_ofst = [&ofst](
      const StateId s,
      const std::set<Label>& valid_labels,
      const std::set<StateId>& valid_states) {
    if (valid_labels.empty() || valid_states.empty()) return;
    auto sit = valid_states.begin();
    std::advance(sit, rand() % valid_states.size());
    auto lit = valid_labels.begin();
    std::advance(lit, rand() % valid_labels.size());
    ofst->AddArc(s, Arc(*lit, *lit, Weight((rand() % 5) / 4.0), *sit));
  };

  auto add_arcs_to_ofst = [&ofst, add_arc_to_ofst](
      const StateId s,
      const int num_arcs,
      const std::set<Label>& valid_labels,
      const std::set<StateId>& valid_states) {
    for (int n = 0; n < num_arcs; ++n) {
      add_arc_to_ofst(s, valid_labels, valid_states);
    }
  };

  auto filter_valid_states = [](
      const StateId s, const std::set<StateId>& all_states) {
    std::set<StateId> valid_states;
    for (StateId v : all_states) {
      if (v > s) { valid_states.emplace(v); }
    }
    return valid_states;
  };

  for (StateId s = 0; s < num_states - 1; ++s) {
    const auto num_arcs = 1 + rand() % 4; // number of output arcs from state s

    auto valid_states_with_delim =
        filter_valid_states(s, states_with_input_delim);
    auto valid_states_with_regular =
        filter_valid_states(s, states_with_input_regular);

    if (states_with_input_delim.count(s) && states_with_input_regular.count(s)) {
      // state s has both types of inputs
      if (ofst->Final(s) != Weight::Zero()) {
        // If the state is final, add all outgoing arcs as delimiters
        add_arcs_to_ofst(s, num_arcs, *delimiters, valid_states_with_delim);
      } else {
        // If the state is not final, we can output one of the two classes,
        // but all outgoing arcs must have the same class.
        if (rand() % 2 == 0 && !valid_states_with_regular.empty()) {
          add_arcs_to_ofst(s, num_arcs, regular, valid_states_with_regular);
        } else {
          add_arcs_to_ofst(s, num_arcs, *delimiters, valid_states_with_delim);
        }
      }
    } else {
      // We can add any type of outgoing label.
      for (int a = 0; a < num_arcs; ++a) {
        if (rand() % 2 == 0 && !valid_states_with_regular.empty()) {
          add_arc_to_ofst(s, regular, valid_states_with_regular);
        } else {
          add_arc_to_ofst(s, *delimiters, valid_states_with_delim);
        }
      }
    }
  }

  Connect(ofst);
}

void TestExpandSubpathsBetweenDelimitersGeneral() {
  typedef StdArc::Label Label;
  typedef StdArc::StateId StateId;
  typedef StdArc::Weight Weight;
  ExpandSubpathsOptions expand_opts;
  VectorFst<StdArc> ifst, ofst;
  // Add four states [0, 1, 2, 3]
  for (size_t n = 0; n < 4; ++n) { ifst.AddState(); }
  // Connect each node to the next, with four labels [0, 1, 2, 3]
  for (StateId s = 0; s < 3; ++s) {
    for (Label label = 0; label <= 3; ++label) {
      ifst.AddArc(s, StdArc(label, label, Weight(s), s + 1));
    }
  }
  ifst.SetStart(0);
  ifst.SetFinal(1, Weight::One());
  ifst.SetFinal(3, Weight::One());

  std::vector<std::tuple<float, std::string, std::string>> expected_paths,
      actual_paths;

  {
    const std::set<int> delimiters = {2};
    auto f = [&delimiters](Label label) -> int {
      return label == 0 ? 0 : (delimiters.count(label) == 0 ? 1 : 2);
    };

    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        f, ifst, &expected_paths, false, delimiters));
    ExpandSubpathsBetweenDelimiters<StdArc>(delimiters, ifst, &ofst);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsWithSameLabelClassBasic[2]", ifst, f,
        false, expected_paths, actual_paths, delimiters);
  }

  {
    const std::set<int> delimiters = {2, 3};
    auto f = [&delimiters](Label label) -> int {
      return label == 0 ? 0 : (delimiters.count(label) == 0 ? 1 : 2);
    };

    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        f, ifst, &expected_paths, false, delimiters));
    ExpandSubpathsBetweenDelimiters<StdArc>(delimiters, ifst, &ofst);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsWithSameLabelClassBasic[2,3]", ifst, f,
        false, expected_paths, actual_paths, delimiters);
  }

  {
    const std::set<int> delimiters = {1, 2, 3};
    auto f = [&delimiters](Label label) -> int {
      return label == 0 ? 0 : (delimiters.count(label) == 0 ? 1 : 2);
    };

    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        f, ifst, &expected_paths, false, std::set<int>{2}));
    ExpandSubpathsBetweenDelimiters<StdArc>(delimiters, ifst, &ofst);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsWithSameLabelClassBasic[1,2,3]", ifst, f,
        false, expected_paths, actual_paths, std::set<int>{2});
  }
}

void TestExpandSubpathsBetweenDelimitersSpecial() {
  typedef StdArc::Label Label;
  typedef StdArc::StateId StateId;
  typedef StdArc::Weight Weight;
  VectorFst<StdArc> ifst, ofst;
  std::vector<std::tuple<float, std::string, std::string>> expected_paths,
      actual_paths;

  for (int r = 0; r < 300; ++r) {
    std::set<Label> delimiters;
    RandFstSpecial(&ifst, &delimiters);

    auto f = [&delimiters](Label label) -> int {
      return label == 0 ? 0 : (delimiters.count(label) == 0 ? 1 : 2);
    };

    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        f, ifst, &expected_paths, false, std::set<int>{2}));
    ExpandSubpathsBetweenDelimiters<StdArc>(delimiters, ifst, &ofst);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsBetweenDelimitersSpecial", ifst, f,
        false, expected_paths, actual_paths, std::set<int>{2});
  }
}

int main() {
  TestExpandSubpathsBetweenDelimitersGeneral();
  TestExpandSubpathsBetweenDelimitersSpecial();
  std::cerr << "Test OK" << std::endl;
}

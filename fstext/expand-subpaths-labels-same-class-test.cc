#undef NDEBUG

#include "fstext/expand-subpaths-labels-same-class.h"
#include "fstext/expand-subpaths-test-utils.h"
#include "fstext/fst-test-utils2.h"
#include "fstext/rand-fst.h"
#include "fstext/kaldi-fst-io.h"

using namespace fst;

void TestExpandSubpathsLabelsSameClassAllSame() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;
  std::set<bool> non_expandable_classes;

  for (int r = 0; r < 200; ++r) {
    RandFstOptions rand_opts;
    rand_opts.acyclic = true;
    rand_opts.n_syms = 2 + rand() % 10;
    rand_opts.n_states = 5 + rand() % 20;
    rand_opts.n_arcs = 10 + kaldi::Rand() % 40;
    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector<std::tuple<float, std::string, std::string>> expected_paths,
        actual_paths;

    auto f = [](Label label) -> bool { return true; };

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<bool>(
        f, *ifst, &expected_paths, true));
    ExpandSubpathsLabelsSameClass<StdArc, bool>(
        f, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, bool>(
        "TestExpandSubpathsLabelsSameClassAllSame.UseInput", *ifst, f,
        true, expected_paths, actual_paths);


    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<bool>(
        f, *ifst, &expected_paths, false));
    ExpandSubpathsLabelsSameClass<StdArc, bool>(
        f, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, bool>(
        "TestExpandSubpathsLabelsSameClassllSame.UseOutput", *ifst, f,
        false, expected_paths, actual_paths);

    delete ifst;
  }
}

void TestExpandSubpathsLabelsSameClassAllDifferent() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;
  std::set<Label> non_expandable_classes;

  for (int r = 0; r < 200; ++r) {
    RandFstOptions rand_opts;
    rand_opts.acyclic = true;
    rand_opts.n_syms = 2 + rand() % 10;
    rand_opts.n_states = 5 + rand() % 20;
    rand_opts.n_arcs = 10 + kaldi::Rand() % 40;
    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector<std::tuple<float, std::string, std::string>>
        expected_paths, actual_paths;

    auto f = [](Label label) -> Label { return label; };

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<Label>(
        f, *ifst, &expected_paths, true));
    ExpandSubpathsLabelsSameClass<StdArc, Label>(
        f, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, Label>(
        "TestExpandSubpathsLabelsSameClassAllDifferent.UseInput", *ifst, f,
        true, expected_paths, actual_paths);

    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<Label>(
        f, *ifst, &expected_paths, false));
    ExpandSubpathsLabelsSameClass<StdArc, Label>(
        f, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, Label>(
        "TestExpandSubpathsLabelsSameClassAllDifferent.UseOutput", *ifst, f,
        false, expected_paths, actual_paths);
    delete ifst;
  }
}

void TestExpandSubpathsLabelsSameClassRandom() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;
  std::set<int> non_expandable_classes;

  for (int r = 0; r < 200; ++r) {
    RandFstOptions rand_opts;
    rand_opts.acyclic = true;
    rand_opts.n_syms = 2 + rand() % 10;
    rand_opts.n_states = 5 + rand() % 20;
    rand_opts.n_arcs = 10 + kaldi::Rand() % 40;
    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector<std::tuple<float, std::string, std::string>>
        expected_paths, actual_paths;
    RandomClass<Label> func;

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        func, *ifst, &expected_paths, true));
    ExpandSubpathsLabelsSameClass<StdArc, int>(
        func, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsLabelsSameClassRandom", *ifst, func, true,
        expected_paths, actual_paths);


    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        func, *ifst, &expected_paths, false));
    ExpandSubpathsLabelsSameClass<StdArc, int>(
        func, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsLabelsSameClassRandom", *ifst, func, false,
        expected_paths, actual_paths);
    delete ifst;
  }
}

void TestExpandSubpathsLabelsSameClassNonExpandableRandom() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;

  for (int r = 0; r < 200; ++r) {
    RandFstOptions rand_opts;
    rand_opts.acyclic = true;
    rand_opts.n_syms = 2 + rand() % 10;
    rand_opts.n_states = 5 + rand() % 20;
    rand_opts.n_arcs = 10 + kaldi::Rand() % 40;

    std::set<int> non_expandable_classes;
    for (int n = 0; n < 3; ++n) {
      non_expandable_classes.emplace(
          (rand() % 2 == 0 ? -1 : +1) * (1 + rand() % 3)
      );
    }

    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector<std::tuple<float, std::string, std::string>>
        expected_paths, actual_paths;
    RandomClass<Label> func;

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        func, *ifst, &expected_paths, true, non_expandable_classes));
    ExpandSubpathsLabelsSameClass<StdArc, int>(
        func, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsLabelsSameClassNonExpandableRandom", *ifst, func,
        true, expected_paths, actual_paths, non_expandable_classes);


    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst<int>(
        func, *ifst, &expected_paths, false, non_expandable_classes));
    ExpandSubpathsLabelsSameClass<StdArc, int>(
        func, *ifst, &ofst, non_expandable_classes, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths<StdArc, int>(
        "TestExpandSubpathsLabelsSameClassNonExpandableRandom", *ifst, func,
        false, expected_paths, actual_paths, non_expandable_classes);
    delete ifst;
  }
}

int main(int argc, char **argv) {
  srand(12345);
  TestExpandSubpathsLabelsSameClassAllSame();
  TestExpandSubpathsLabelsSameClassAllDifferent();
  TestExpandSubpathsLabelsSameClassRandom();
  TestExpandSubpathsLabelsSameClassNonExpandableRandom();
  std::cerr << "Test OK" << std::endl;
}

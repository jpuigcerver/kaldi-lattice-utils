#undef NDEBUG

#include "fst/script/compile.h"
#include "fst/script/print.h"

#include "fstext/expand-subpaths-same-label-class.h"
#include "fstext/fstext-utils.h"
#include "fstext/fst-test-utils2.h"
#include "fstext/rand-fst.h"
#include "fstext/kaldi-fst-io.h"

using namespace fst;

template<typename F, typename Label>
void SplitPathByLabelClass(
    const F& f,
    const std::vector<Label> &isyms,
    const std::vector<Label> &osyms,
    std::vector<std::string> *isubpaths,
    std::vector<std::string> *osubpaths,
    bool use_input) {
  KALDI_ASSERT(isyms.size() == osyms.size());
  isubpaths->clear();
  osubpaths->clear();
  const auto cls_eps = f(0);

  auto cls_prev = cls_eps;
  std::vector<Label> isubpath, osubpath;
  for (size_t i = 0; i < isyms.size(); ++i) {
    const auto cls_curr = f(use_input ? isyms[i] : osyms[i]);
    if (cls_curr == cls_prev || cls_prev == cls_eps || cls_curr == cls_eps) {
      isubpath.push_back(isyms[i]);
      osubpath.push_back(osyms[i]);
    } else {
      const auto istr = LabelVectorToString(isubpath, nullptr, kStringSeparator);
      const auto ostr = LabelVectorToString(osubpath, nullptr, kStringSeparator);
      if (!istr.empty()) isubpaths->emplace_back(istr);
      if (!ostr.empty()) osubpaths->emplace_back(ostr);
      isubpath.clear();
      osubpath.clear();
      isubpath.push_back(isyms[i]);
      osubpath.push_back(osyms[i]);
    }
    cls_prev = cls_curr != cls_eps ? cls_curr : cls_prev;
  }

  const auto istr = LabelVectorToString(isubpath, nullptr, kStringSeparator);
  const auto ostr = LabelVectorToString(osubpath, nullptr, kStringSeparator);
  if (!istr.empty()) isubpaths->emplace_back(istr);
  if (!ostr.empty()) osubpaths->emplace_back(ostr);
}

std::string StringVectorToString(const std::vector<std::string> &v,
                                 const std::string &sep = " ") {
  std::string str;
  for (const auto &s : v) {
    if (!str.empty()) str += sep;
    str += s;
  }
  return str;
}

template <typename F>
bool ExpandAndGetPathsOfSubpathsFromFst(
    const F& f,
    const VectorFst <StdArc> &vfst,
    std::vector<std::tuple<float, std::string, std::string>> *paths,
    bool use_input) {
  typedef StdArc::Label Label;

  if (vfst.Properties(kAcyclic, true) != kAcyclic)
    return false;

  std::vector<VectorFst<StdArc>> vtmp;
  {
    VectorFst <StdArc> tmp;
    ShortestPath(vfst, &tmp, std::numeric_limits<int>::max());
    ConvertNbestToVector(tmp, &vtmp);
  }

  paths->clear();
  for (const VectorFst <StdArc> &path : vtmp) {
    std::vector<Label> path_isyms, path_osyms;
    StdArc::Weight path_cost;
    GetLinearSymbolSequence(path, &path_isyms, &path_osyms, &path_cost, /*include_eps=*/true);

    std::vector<std::string> isubpaths, osubpaths;
    SplitPathByLabelClass(f, path_isyms, path_osyms, &isubpaths, &osubpaths, use_input);

    paths->emplace_back(path_cost.Value(), StringVectorToString(isubpaths), StringVectorToString(osubpaths));
  }
  std::sort(paths->begin(), paths->end());

  return true;
}

bool GetPathsOfSubpathsFromExpandedFst(
    const VectorFst <StdArc> &vfst,
    std::vector <std::tuple<float, std::string, std::string>> *paths) {
  if (vfst.Properties(kAcyclic, true) != kAcyclic)
    return false;

  std::vector <VectorFst<StdArc>> vtmp;
  {
    VectorFst <StdArc> tmp;
    ShortestPath(vfst, &tmp, std::numeric_limits<int>::max());
    ConvertNbestToVector(tmp, &vtmp);
  }

  paths->clear();
  for (const VectorFst <StdArc> &path : vtmp) {
    std::vector <StdArc::Label> path_isyms, path_osyms;
    StdArc::Weight path_cost;
    GetLinearSymbolSequence(path, &path_isyms, &path_osyms, &path_cost);

    const auto
        isubpaths = LabelVectorToString(path_isyms, vfst.InputSymbols());
    const auto
        osubpaths = LabelVectorToString(path_osyms, vfst.OutputSymbols());
    paths->emplace_back(path_cost.Value(), isubpaths, osubpaths);
  }
  std::sort(paths->begin(), paths->end());

  return true;
}

std::string PathsOfSubpathsToString(
    const std::vector <std::tuple<float, std::string, std::string>> &paths) {
  std::string str;
  for (size_t i = 0; i < paths.size(); ++i) {
    const auto &p = paths[i];
    str += std::to_string(i + 1) + " / " + std::to_string(std::get<0>(p))
        + " / " + std::get<1>(p) + " / " + std::get<2>(p) + "\n";
  }
  return str;
}

template <typename Arc, typename F>
std::string LabelClassFunctionSummary(const Fst<Arc>& ifst, const F& func, bool use_input) {
  std::set<typename Arc::Label> labels{0};
  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      labels.insert(use_input ? ait.Value().ilabel : ait.Value().olabel);
    }
  }

  std::ostringstream oss;
  auto it = labels.begin();
  oss << *it << ": " << func(*it);
  for (++it; it != labels.end(); ++it) {
    oss << ", " << *it << ": " << func(*it);
  }
  return oss.str();
}


template<typename Arc, typename F>
void TestPathsOfSubpaths(
    const std::string &test_name,
    const Fst<Arc> &ifst,
    const F& func,
    const bool use_input,
    const std::vector <std::tuple<float, std::string, std::string>> &expected,
    const std::vector <std::tuple<float, std::string, std::string>> &actual) {
  if (!AreBestPathsEqual(expected, actual)) {
    std::ostringstream oss;
    fst::script::PrintFst(
        ifst,
        oss,
        "",
        ifst.InputSymbols(),
        ifst.OutputSymbols());
    KALDI_ERR << "FAILED TEST " + test_name + ":\n"
              << "Original FST:\n"
              << oss.str() << "\n"
              << "Label classes:\n"
              << LabelClassFunctionSummary(ifst, func, use_input) << "\n"
              << "Expected paths:\n"
              << PathsOfSubpathsToString(expected) << "\n"
              << "Actual paths:\n"
              << PathsOfSubpathsToString(actual) << "\n";
  }
}

void TestExpandSubpathsWithSameLabelClassAllSame() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;
  RandFstOptions rand_opts;
  rand_opts.acyclic = true;

  for (int r = 0; r < 100; ++r) {
    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector <std::tuple<float, std::string, std::string>> expected_paths, actual_paths;

    auto f = [](Label label) -> bool { return true; };

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst(f, *ifst, &expected_paths, true));
    const auto num_arcs1 = ExpandSubpathsWithSameLabelClass<StdArc, bool>(f, *ifst, &ofst, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths("TestExpandSubpathsWithSameLabelClassAllSame", *ifst, f, true, expected_paths, actual_paths);
    KALDI_ASSERT(num_arcs1 == expected_paths.size());


    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst(f, *ifst, &expected_paths, false));
    const auto num_arcs2 = ExpandSubpathsWithSameLabelClass<StdArc, bool>(f, *ifst, &ofst, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths("TestExpandSubpathsWithSameLabelClassAllSame", *ifst, f, false, expected_paths, actual_paths);
    KALDI_ASSERT(num_arcs2 == expected_paths.size());
    delete ifst;
  }
}

void TestExpandSubpathsWithSameLabelClassAllDifferent() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;
  RandFstOptions rand_opts;
  rand_opts.acyclic = true;

  for (int r = 0; r < 100; ++r) {
    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector <std::tuple<float, std::string, std::string>> expected_paths, actual_paths;

    auto f = [](Label label) -> Label { return label; };

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst(f, *ifst, &expected_paths, true));
    ExpandSubpathsWithSameLabelClass<StdArc, Label>(f, *ifst, &ofst, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths("TestExpandSubpathsWithSameLabelClassAllDifferent", *ifst, f, true, expected_paths, actual_paths);


    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst(f, *ifst, &expected_paths, false));
    ExpandSubpathsWithSameLabelClass<StdArc, Label>(f, *ifst, &ofst, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths("TestExpandSubpathsWithSameLabelClassAllDifferent", *ifst, f, false, expected_paths, actual_paths);
    delete ifst;
  }
}


void TestExpandSubpathsWithSameLabelClassRandom() {
  typedef StdArc::Label Label;
  ExpandSubpathsOptions expand_opts;
  RandFstOptions rand_opts;
  rand_opts.acyclic = true;

  for (int r = 0; r < 100; ++r) {
    VectorFst<StdArc> *ifst = RandFst<StdArc>(rand_opts), ofst;
    std::vector <std::tuple<float, std::string, std::string>> expected_paths, actual_paths;
    RandomClass<Label> func;

    // use_input = true
    expand_opts.use_input = true;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst(func, *ifst, &expected_paths, true));
    ExpandSubpathsWithSameLabelClass<StdArc, int>(func, *ifst, &ofst, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths("TestExpandSubpathsWithSameLabelClassRandom", *ifst, func, true, expected_paths, actual_paths);


    // use_input = false
    expand_opts.use_input = false;
    KALDI_ASSERT(ExpandAndGetPathsOfSubpathsFromFst(func, *ifst, &expected_paths, false));
    ExpandSubpathsWithSameLabelClass<StdArc, int>(func, *ifst, &ofst, expand_opts);
    KALDI_ASSERT(GetPathsOfSubpathsFromExpandedFst(ofst, &actual_paths));
    TestPathsOfSubpaths("TestExpandSubpathsWithSameLabelClassRandom", *ifst, func, false, expected_paths, actual_paths);
    delete ifst;
  }
}

int main(int argc, char** argv) {
  InitFst("USAGE", &argc, &argv, true);
  srand(12345);
  TestExpandSubpathsWithSameLabelClassAllSame();
  TestExpandSubpathsWithSameLabelClassAllDifferent();
  TestExpandSubpathsWithSameLabelClassRandom();
  std::cerr << "Test OK" << std::endl;
}

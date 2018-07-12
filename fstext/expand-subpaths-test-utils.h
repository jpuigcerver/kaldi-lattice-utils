#ifndef KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_TEST_UTILS_H_
#define KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_TEST_UTILS_H_

#include "fst/script/print.h"
#include "fstext/fstext-utils.h"
#include "fstext/fst-test-utils2.h"

namespace fst {

template<typename ClassType, typename F, typename Label>
void SplitPathByLabelClass(
    const F &f,
    const std::vector<Label> &isyms,
    const std::vector<Label> &osyms,
    std::vector<std::string> *isubpaths,
    std::vector<std::string> *osubpaths,
    bool use_input,
    const std::set<ClassType>& non_expandable_symbols) {
  KALDI_ASSERT(isyms.size() == osyms.size());
  isubpaths->clear();
  osubpaths->clear();
  const auto cls_eps = f(0);

  auto cls_prev = cls_eps;
  std::vector<Label> isubpath, osubpath;
  for (size_t i = 0; i < isyms.size(); ++i) {
    auto cls_curr = f(use_input ? isyms[i] : osyms[i]);
    if (cls_curr == cls_eps) { cls_curr = cls_prev; }
    if (cls_prev == cls_eps ||
        (cls_curr == cls_prev && non_expandable_symbols.count(cls_curr) == 0)) {
      isubpath.push_back(isyms[i]);
      osubpath.push_back(osyms[i]);
    } else {
      const auto istr =
          LabelVectorToString(isubpath, nullptr, kStringSeparator);
      const auto ostr =
          LabelVectorToString(osubpath, nullptr, kStringSeparator);
      if (!istr.empty()) isubpaths->emplace_back(istr);
      if (!ostr.empty()) osubpaths->emplace_back(ostr);
      isubpath.clear();
      osubpath.clear();
      isubpath.push_back(isyms[i]);
      osubpath.push_back(osyms[i]);
    }
    cls_prev = cls_curr;
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

template<typename ClassType, typename F>
bool ExpandAndGetPathsOfSubpathsFromFst(
    const F &f,
    const VectorFst<StdArc> &vfst,
    std::vector<std::tuple<float, std::string, std::string>> *paths,
    bool use_input,
    const std::set<ClassType>& non_expandable_labels = std::set<ClassType>()) {
  typedef StdArc::Label Label;

  if (vfst.Properties(kAcyclic, true) != kAcyclic)
    return false;

  std::vector<VectorFst<StdArc>> vtmp;
  {
    VectorFst<StdArc> tmp;
    ShortestPath(vfst, &tmp, std::numeric_limits<int>::max());
    ConvertNbestToVector(tmp, &vtmp);
  }

  paths->clear();
  for (const VectorFst<StdArc> &path : vtmp) {
    std::vector<Label> path_isyms, path_osyms;
    StdArc::Weight path_cost;
    GetLinearSymbolSequence(
        path, &path_isyms, &path_osyms, &path_cost, /*include_eps=*/true);

    std::vector<std::string> isubpaths, osubpaths;
    SplitPathByLabelClass<ClassType, F, StdArc::Label>(
        f, path_isyms, path_osyms, &isubpaths, &osubpaths, use_input,
        non_expandable_labels);

    paths->emplace_back(
        path_cost.Value(),
        StringVectorToString(isubpaths),
        StringVectorToString(osubpaths));
  }
  std::sort(paths->begin(), paths->end());

  return true;
}

bool GetPathsOfSubpathsFromExpandedFst(
    const VectorFst<StdArc> &vfst,
    std::vector<std::tuple<float, std::string, std::string>> *paths) {
  if (vfst.Properties(kAcyclic, true) != kAcyclic)
    return false;

  std::vector<VectorFst<StdArc>> vtmp;
  {
    VectorFst<StdArc> tmp;
    ShortestPath(vfst, &tmp, std::numeric_limits<int>::max());
    ConvertNbestToVector(tmp, &vtmp);
  }

  paths->clear();
  for (const VectorFst<StdArc> &path : vtmp) {
    std::vector<StdArc::Label> path_isyms, path_osyms;
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
    const std::vector<std::tuple<float, std::string, std::string>> &paths) {
  std::string str;
  for (size_t i = 0; i < paths.size(); ++i) {
    const auto &p = paths[i];
    str += std::to_string(i + 1) + " / " + std::to_string(std::get<0>(p))
        + " / " + std::get<1>(p) + " / " + std::get<2>(p) + "\n";
  }
  return str;
}

template<typename Arc, typename F>
std::string LabelClassFunctionSummary(
    const Fst<Arc> &ifst, const F &func, bool use_input) {
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

template <typename T>
std::string SetToString(const std::set<T>& s) {
  std::ostringstream oss;
  for (const auto& v: s) {
    oss << v << " ";
  }
  return oss.str();
}

template<typename Arc, typename ClassType, typename F>
void TestPathsOfSubpaths(
    const std::string &test_name,
    const Fst<Arc> &ifst,
    const F &func,
    const bool use_input,
    const std::vector<std::tuple<float, std::string, std::string>> &expected,
    const std::vector<std::tuple<float, std::string, std::string>> &actual,
    const std::set<ClassType>& non_expandable_classes = std::set<ClassType>()) {
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
              << "Non-expandable classes:\n"
              << SetToString(non_expandable_classes) << "\n"
              << "Expected paths:\n"
              << PathsOfSubpathsToString(expected) << "\n"
              << "Actual paths:\n"
              << PathsOfSubpathsToString(actual) << "\n";
  } else {
    /*
    // USED ONLY FOR DEBUGGING
    std::ostringstream oss;
    fst::script::PrintFst(
        ifst,
        oss,
        "",
        ifst.InputSymbols(),
        ifst.OutputSymbols());
    KALDI_LOG << "SUCCEEDED TEST " + test_name + ":\n"
              << "Original FST:\n"
              << oss.str() << "\n"
              << "Label classes:\n"
              << LabelClassFunctionSummary(ifst, func, use_input) << "\n"
              << "Non-expandable classes:\n"
              << SetToString(non_expandable_classes) << "\n"
              << "Expected paths:\n"
              << PathsOfSubpathsToString(expected) << "\n"
              << "Actual paths:\n"
              << PathsOfSubpathsToString(actual) << "\n";
    */
  }
}

}  // namespace fst

#endif //KALDI_LATTICE_UTILS_FSTEXT_EXPAND_SUBPATHS_TEST_UTILS_H_

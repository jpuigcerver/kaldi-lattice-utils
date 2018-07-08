#undef NDEBUG

#include "fstext/make-preceding-symbols-same-class.h"
#include "fstext/rand-fst.h"
#include "fst/script/print.h"
#include "fstext/fstext-utils.h"
#include "fstext/fst-test-utils2.h"

using namespace fst;

template <typename StateId, typename ClassType>
inline bool CheckStateClass(
    const StateId& state,
    const ClassType& cls,
    const std::vector<ClassType>& state_class) {
  if (state >= state_class.size()) {
    KALDI_WARN << "State " << state << " not present in state_class (size = "
               << state_class.size() << ")";
    return false;
  }
  if (state_class[state] != cls) {
    KALDI_WARN << "State " << state << " class does not match "
               << "(expected = " << state_class[state]
               << ", found = " << cls << ")";
    return false;
  }
  return true;
}

template <typename Arc, typename ClassType, typename F>
bool CheckPrecedingSymbolsAreSameClass(
    const Fst<Arc>& ifst,
    const F &f,
    const std::vector<ClassType>& state_class,
    const PrecedingSymbolsSameClassOptions& opts = PrecedingSymbolsSameClassOptions()) {
  if (ifst.Start() == kNoStateId) { return state_class.empty(); }

  const auto c_eps = f(0);

  if (!CheckStateClass(ifst.Start(), c_eps, state_class))
    return false;

  for (StateIterator<Fst<Arc>> sit(ifst); !sit.Done(); sit.Next()) {
    const auto s = sit.Value();
    const auto c = state_class[s];
    for (ArcIterator<Fst<Arc>> ait(ifst, s); !ait.Done(); ait.Next()) {
      const auto& arc = ait.Value();
      const auto c_arc = f(opts.use_input ? arc.ilabel : arc.olabel);
      const auto c_nstate =
          opts.propagate_epsilon_class ? (c_arc != c_eps ? c_arc : c) : c_arc;
      if (!CheckStateClass(arc.nextstate, c_nstate, state_class)) {
        return false;
      }
    }
  }

  return true;
}

template <typename Arc>
void TestMakePrecedingSymbolsSameClassAllEqual() {
  // Test MakePrecedingSymbolsSameClass when all input labels belong to the
  // same class. In this case, the output FST is equal to the input FST,
  // except the numbering of the states.
  // Note: propagate_epsilon_class does not have any effect here
  typedef typename Arc::Label Label;
  for (int r = 0; r < 100; ++r) {
    VectorFst<Arc> *ifst = RandFst<Arc>(), ofst;
    ArcSort(ifst, ILabelCompare<Arc>());
    RenumberStatesFst(ifst);

    std::vector<bool> state_class;
    PrecedingSymbolsSameClassOptions opts;

    // use_input = false, propagate_epsilon_class = false
    MakePrecedingSymbolsSameClass([](Label label) { return true; }, *ifst, &ofst, &state_class, opts);
    RenumberStatesFst(&ofst);
    KALDI_ASSERT(Equal(*ifst, ofst));
    KALDI_ASSERT(state_class == std::vector<bool>(ofst.NumStates(), true));

    // use_input = true, propagate_epsilon_class = false
    opts.use_input = true;
    MakePrecedingSymbolsSameClass([](Label label) { return true; }, *ifst, &ofst, &state_class, opts);
    RenumberStatesFst(&ofst);
    KALDI_ASSERT(Equal(*ifst, ofst));
    KALDI_ASSERT(state_class == std::vector<bool>(ofst.NumStates(), true));

    delete ifst;
  }
}

template <typename Arc>
void TestMakePrecedingSymbolsSameClassAllDifferent() {
  // Test MakePrecedingSymbolsSameClass when all input labels belong to
  // different classes.
  typedef typename Arc::Label Label;
  RandFstOptions rand_opts;
  rand_opts.acyclic = true;  // We need acyclic = false to compare ALL paths

  for (int r = 0; r < 100; ++r) {
    VectorFst<Arc> *ifst = RandFst<Arc>(rand_opts), ofst;

    std::vector<Label> state_class;
    PrecedingSymbolsSameClassOptions opts;

    // use_input = false, propagate_epsilon_class = false
    opts.use_input = false; opts.propagate_epsilon_class = false;
    MakePrecedingSymbolsSameClass(
        [](Label label) { return label; }, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());
    KALDI_ASSERT(CheckPrecedingSymbolsAreSameClass(
        ofst, [](Label label) { return label; }, state_class, opts));


    // use_input = true, propagate_epsilon_class = false
    opts.use_input = true; opts.propagate_epsilon_class = false;
    MakePrecedingSymbolsSameClass(
        [](Label label) { return label; }, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());
    KALDI_ASSERT(CheckPrecedingSymbolsAreSameClass(
        ofst, [](Label label) { return label; }, state_class, opts));

    // use_input = false, propagate_epsilon_class = true
    opts.use_input = false; opts.propagate_epsilon_class = true;
    MakePrecedingSymbolsSameClass(
        [](Label label) { return label; }, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());
    KALDI_ASSERT(CheckPrecedingSymbolsAreSameClass(
        ofst, [](Label label) { return label; }, state_class, opts));

    // use_input = true, propagate_epsilon_class = true
    opts.use_input = true; opts.propagate_epsilon_class = true;
    MakePrecedingSymbolsSameClass(
        [](Label label) { return label; }, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());
    KALDI_ASSERT(CheckPrecedingSymbolsAreSameClass(
        ofst, [](Label label) { return label; }, state_class, opts));
  }
}

template <typename Arc>
void TestMakePrecedingSymbolsSameClassRandom() {
  typedef typename Arc::Label Label;
  RandFstOptions rand_opts;
  rand_opts.acyclic = true;  // We need acyclic = false to compare ALL paths

  for (int r = 0; r < 100; ++r) {
    VectorFst<Arc> *ifst = RandFst<Arc>(rand_opts), ofst;
    RandomClass<Label> func;

    std::vector<int> state_class;
    PrecedingSymbolsSameClassOptions opts;

    // use_input = false, propagate_epsilon_class = false
    opts.use_input = false; opts.propagate_epsilon_class = false;
    MakePrecedingSymbolsSameClass(func, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());

    // use_input = false, propagate_epsilon_class = false
    opts.use_input = false; opts.propagate_epsilon_class = true;
    MakePrecedingSymbolsSameClass(func, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());

    // use_input = true, propagate_epsilon_class = false
    opts.use_input = true; opts.propagate_epsilon_class = false;
    MakePrecedingSymbolsSameClass(func, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());

    // use_input = true, propagate_epsilon_class = true
    opts.use_input = true; opts.propagate_epsilon_class = true;
    MakePrecedingSymbolsSameClass(func, *ifst, &ofst, &state_class, opts);
    KALDI_ASSERT(AreBestPathsEqual(*ifst, ofst));
    KALDI_ASSERT(state_class.size() == ofst.NumStates());
  }
}

int main() {
  srand(12345);
  TestMakePrecedingSymbolsSameClassAllEqual<StdArc>();
  TestMakePrecedingSymbolsSameClassAllDifferent<StdArc>();
  TestMakePrecedingSymbolsSameClassRandom<StdArc>();
  std::cerr << "Test OK" << std::endl;
  return 0;
}

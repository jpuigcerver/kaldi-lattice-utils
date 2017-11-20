#undef NDEBUG   // Undefine NDEBUG, or tests may fail.

#include "util/kaldi-io.h"
#include "fstext/fstext-utils2.h"

using namespace kaldi;

void BuildFst(fst::VectorFst<fst::StdArc>* ofst) {
  typedef fst::StdArc::Weight Weight;
  ofst->DeleteStates();

  ofst->AddState();
  ofst->AddState();
  ofst->AddState();

  ofst->SetStart(0);
  ofst->AddArc(0, fst::StdArc(0, 0, Weight::One(), 1));
  ofst->AddArc(0, fst::StdArc(1, 2, Weight::One(), 1));
  ofst->AddArc(1, fst::StdArc(2, 1, Weight::One(), 2));
  ofst->AddArc(1, fst::StdArc(1, 2, Weight::One(), 2));
  ofst->AddArc(1, fst::StdArc(4, 3, Weight::One(), 2));
  ofst->AddArc(2, fst::StdArc(3, 4, Weight::One(), 2));
  ofst->SetFinal(2, Weight::One());
}

void DisambiguateStatesByInputLabelGroup() {
  typedef fst::StdArc::Label Label;
  fst::VectorFst<fst::StdArc> tfst;
  BuildFst(&tfst);

  // Note: Labels 3 and 4 are assigned to a new group (ID = MAX_INT).
  const std::unordered_map<Label, int> label_group_2 =
      {{0, 0}, {1, 1}, {2, 1}};
  {
    // Note: All labels are assigned to separate groups.
    const std::unordered_map<Label, int> label_group =
      {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};
    fst::VectorFst<fst::StdArc> ofst;
    std::vector<int> state_group;
    DisambiguateStatesByInputLabelGroup(tfst, label_group, &ofst, &state_group,
                                        true);
  }
}

int main(int argc, char** argv) {
  DisambiguateStatesByInputLabelGroup();
  std::cout << "Test OK." << std::endl;
}

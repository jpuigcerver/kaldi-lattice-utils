#undef NDEBUG

#include "fst/fstlib.h"
#include "fstext/label-group.h"

using namespace fst;

void TestAddGroup() {
  LabelGroup<int> label_group;
  KALDI_ASSERT(label_group.NumGroups() == 1);  // Epsilon's groups
  label_group.AddGroup({1, 3, 5});
  label_group.AddGroup({2, 4, 6});
  KALDI_ASSERT(label_group.NumGroups() == 3);
}

void TestGetLabelGroup() {
  LabelGroup<int> label_group;
  KALDI_ASSERT(label_group[0] == 0);

  label_group.AddGroup({1, 3, 5});
  label_group.AddGroup({2, 4, 6});
  KALDI_ASSERT(label_group[0] == 0);
  KALDI_ASSERT(label_group[1] == 1);
  KALDI_ASSERT(label_group[3] == 1);
  KALDI_ASSERT(label_group[5] == 1);
  KALDI_ASSERT(label_group[2] == 2);
  KALDI_ASSERT(label_group[4] == 2);
  KALDI_ASSERT(label_group[6] == 2);
  KALDI_ASSERT(label_group[99] == 3);

  for (int label = 0; label < 10; ++label) {
    KALDI_ASSERT(label_group[label] == label_group(label));
  }
}

void TestCopyConstructor() {
  LabelGroup<int> label_group1;
  label_group1.AddGroup({1, 3, 5});
  LabelGroup<int> label_group2(label_group1);
  KALDI_ASSERT(label_group1.NumGroups() == label_group2.NumGroups());
  for (int label = 0; label < 10; ++label) {
    KALDI_ASSERT(label_group1(label) == label_group2(label));
  }
}

void TestParseString() {
  std::vector<std::string> groups_str{
    "1 3 5 ; 2 4 6",
    "; 1 3 5 ; 2 4 6 ;",
    ";1 3 5;2 4 6;",
    ";;;;   1   3    5;   2 4     6;",
  };
  for (const auto& str : groups_str) {
    LabelGroup<int> label_group;
    label_group.ParseString(str);
    KALDI_ASSERT(label_group.NumGroups() == 3);
    KALDI_ASSERT(label_group(0) == 0);
    KALDI_ASSERT(label_group(1) == 1);
    KALDI_ASSERT(label_group(3) == 1);
    KALDI_ASSERT(label_group(5) == 1);
    KALDI_ASSERT(label_group(2) == 2);
    KALDI_ASSERT(label_group(4) == 2);
    KALDI_ASSERT(label_group(6) == 2);
    KALDI_ASSERT(label_group(7) == 3);
  }
}

int main() {
  TestAddGroup();
  TestGetLabelGroup();
  TestCopyConstructor();
  TestParseString();
  std::cerr << "Test OK" << std::endl;
  return 0;
}

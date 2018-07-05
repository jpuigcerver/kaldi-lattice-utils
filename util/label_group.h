#ifndef KALDI_LATTICE_UTILS_LABEL_GROUP_H_
#define KALDI_LATTICE_UTILS_LABEL_GROUP_H_

#include <unordered_map>
#include <vector>

#include "util/text-utils.h"

namespace kaldi {

template <typename Label>
class LabelGroup {
 public:
  LabelGroup() {
    map_[0] = 0;
    num_groups_ = 1;
  }

  void AddGroup(const Container& labels) {
    for (const auto& label : labels) {
      auto r = label_group->emplace(label, num_groups_);
      if (!r.second) {
        KALDI_ERR << "Each label must be assigned to one group at most. "
                  << "Label " << label << " was assigned to both groups "
                  << r.first->second << " and " << num_groups_ << ".";
      }
    }
    ++num_groups_;
  }

  void AddGroupsStr(const std::string& multi_group_str) {
    std::vector<std::string> tmp;
    kaldi::SplitStringToVector(separator_groups_str, ";", true, &tmp);
    for (const auto& group_str : tmp) {
      std::vector<Label> tmp;
      kaldi::SplitStringToIntegers(tmp[i], " ", true, &tmp2);
      AddGroup(tmp2);
    }
  }

  Label operator[](const Label& label) const {
    auto it = map_.find(label);
    if (it == map_.end()) {
      return std::numeric_limits<Label>::max();
    } else {
      return it->second;
    }
  }

  Label NumGroups() const {
    return num_groups_;
  }

 protected:
  std::unordered_map<Label, Label> map_;
  Label num_groups_;
};

}  // namespace kaldi

#endif //KALDI_LATTICE_UTILS_LABEL_GROUP_H_

#ifndef KALDI_LATTICE_UTILS_LABEL_GROUP_H_
#define KALDI_LATTICE_UTILS_LABEL_GROUP_H_

#include <unordered_map>
#include <vector>

#include "util/text-utils.h"

namespace fst {

template <typename Label>
class LabelGroup {
 public:
  LabelGroup() {
    map_[0] = 0;
    num_groups_ = 1;
  }

  LabelGroup(const LabelGroup& other)
      : map_(other.map_), num_groups_(other.num_groups_) {}

  void AddGroup(const std::vector<Label>& labels) {
    for (const auto& label : labels) {
      auto r = map_.emplace(label, num_groups_);
      if (!r.second) {
        KALDI_WARN << "Each label must be assigned to one group at most. "
                   << "Label " << label << " changed from group "
                   << r.first->second << " to " << num_groups_ << ".";
        r.first->second = num_groups_;
      }
    }
    ++num_groups_;
  }

  void ParseString(const std::string& multi_group_str) {
    std::vector<std::string> tmp;
    kaldi::SplitStringToVector(multi_group_str, ";", true, &tmp);
    for (const auto& group_str : tmp) {
      std::vector<Label> tmp2;
      kaldi::SplitStringToIntegers(group_str, " ", true, &tmp2);
      AddGroup(tmp2);
    }
  }

  Label operator[](const Label& label) const {
    auto it = map_.find(label);
    if (it == map_.end()) {
      return num_groups_;
    } else {
      return it->second;
    }
  }

  Label operator()(const Label& label) const {
    return operator[](label);
  }

  Label NumGroups() const {
    return num_groups_;
  }

 protected:
  std::unordered_map<Label, Label> map_;
  Label num_groups_;
};

}  // namespace fst

#endif //KALDI_LATTICE_UTILS_LABEL_GROUP_H_

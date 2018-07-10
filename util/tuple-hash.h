#ifndef KALDI_LATTICE_UTILS_UTIL_TUPLE_HASH_H_
#define KALDI_LATTICE_UTILS_UTIL_TUPLE_HASH_H_

#include <tuple>

namespace kaldi {

template<typename T>
struct hash {
  inline size_t operator()(const T &v) const {
    return std::hash<T>()(v);
  }
};

namespace {

template<class T>
inline void hash_combine(std::size_t &seed, const T &v) {
  seed ^= kaldi::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Recursive template code derived from Matthieu M.
template<class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
struct TupleHashImpl {
  static void apply(size_t &seed, Tuple const &tuple) {
    TupleHashImpl<Tuple, Index - 1>::apply(seed, tuple);
    hash_combine(seed, std::get<Index>(tuple));
  }
};

template<class Tuple>
struct TupleHashImpl<Tuple, 0> {
  static void apply(size_t &seed, Tuple const &tuple) {
    hash_combine(seed, std::get<0>(tuple));
  }
};

}  // namespace


// Specialization of the hash struct for std::tuple
template<typename ... T>
struct hash<std::tuple < T...> > {
inline size_t operator()(const std::tuple<T...> &v) const {
  size_t seed = 0;
  TupleHashImpl<std::tuple < T...> > ::apply(seed, v);
  return seed;
}
};

}  // namespace kaldi

#endif //KALDI_LATTICE_UTILS_UTIL_TUPLE_HASH_H_

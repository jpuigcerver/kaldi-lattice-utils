// fstext/log-tropical-weight.h

// Copyright (c) 2017 Joan Puigcerver <joapuipe@upv.es>

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#ifndef KALDI_FSTEXT_LOG_TROPICAL_WEIGHT_H_
#define KALDI_FSTEXT_LOG_TROPICAL_WEIGHT_H_

namespace fst {

///  LogTropical semiring: (min, -log(e^-x + e^-y), inf, inf).
template <class T>
class LogTropicalWeightTpl : public FloatWeightTpl<T> {
 public:
  using typename FloatWeightTpl<T>::ValueType;
  using FloatWeightTpl<T>::Value;
  using ReverseWeight = LogTropicalWeightTpl<T>;

  constexpr LogTropicalWeightTpl() : FloatWeightTpl<T>() {}

  constexpr LogTropicalWeightTpl(T f) : FloatWeightTpl<T>(f) {}

  constexpr LogTropicalWeightTpl(const LogTropicalWeightTpl<T> &weight)
      : FloatWeightTpl<T>(weight) {}

  static const LogTropicalWeightTpl<T> &Zero() {
    static const LogTropicalWeightTpl zero(FloatLimits<T>::PosInfinity());
    return zero;
  }

  static const LogTropicalWeightTpl<T> &One() {
    static const LogTropicalWeightTpl one(FloatLimits<T>::PosInfinity());
    return one;
  }

  static const LogTropicalWeightTpl<T> &NoWeight() {
    static const LogTropicalWeightTpl no_weight(FloatLimits<T>::NumberBad());
    return no_weight;
  }

  static const string &Type() {
    static const string type =
        string("logtropical") + FloatWeightTpl<T>::GetPrecisionString();
    return type;
  }

  bool Member() const {
    ///  First part fails for IEEE NaN.
    return Value() == Value() && Value() != FloatLimits<T>::NegInfinity();
  }

  LogTropicalWeightTpl<T> Quantize(float delta = kDelta) const {
    if (Value() == FloatLimits<T>::NegInfinity() ||
        Value() == FloatLimits<T>::PosInfinity() || Value() != Value()) {
      return *this;
    } else {
      return LogTropicalWeightTpl<T>(floor(Value() / delta + 0.5F) * delta);
    }
  }

  LogTropicalWeightTpl<T> Reverse() const { return *this; }

  static constexpr uint64 Properties() {
    return kLeftSemiring | kRightSemiring | kCommutative | kPath | kIdempotent;
  }
};

// a + (inf) = a
template <class T>
inline LogTropicalWeightTpl<T> Plus(const LogTropicalWeightTpl<T> &w1,
                                    const LogTropicalWeightTpl<T> &w2) {
  if (!w1.Member() || !w2.Member()) return LogTropicalWeightTpl<T>::NoWeight();
  return w1.Value() < w2.Value() ? w1 : w2;
}

// a * (inf) = (inf)
template <class T>
inline LogTropicalWeightTpl<T> Times(const LogTropicalWeightTpl<T> &w1,
                                     const LogTropicalWeightTpl<T> &w2) {
  T f1 = w1.Value(), f2 = w2.Value();
  if (f1 == FloatLimits<T>::PosInfinity()) {
    return w2;
  } else if (f2 == FloatLimits<T>::PosInfinity()) {
    return w1;
  } else if (f1 > f2) {
    return LogTropicalWeightTpl<T>(f2 - internal::LogPosExp(f1 - f2));
  } else {
     return LogTropicalWeightTpl<T>(f1 - internal::LogPosExp(f2 - f1));
  }
}

}  // namespace fst




#endif  // KALDI_FSTEXT_LOG_TROPICAL_WEIGHT_H_

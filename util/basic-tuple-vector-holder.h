// util/basic-tuple-vector-holder.h

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

#include <vector>
#include <tuple>

#include "util/kaldi-io.h"
#include "util/text-utils.h"

namespace kaldi {

namespace internal {

// Utility classes to write an arbitrary tuple to an output stream

template <class T>
void ReadTupleElement(std::istream& is, bool binary, T* elem) {
  ReadBasicType(is, binary, elem);
}

template <class T>
void WriteTupleElement(std::ostream& os, bool binary, const T& elem) {
  WriteBasicType(os, binary, elem);
}

// Specialization for std::string
template <>
void ReadTupleElement(std::istream& is, bool binary, std::string* elem) {
  ReadToken(is, binary, elem);
}

template <>
void WriteTupleElement<std::string>(
    std::ostream& os, bool binary, const std::string& elem) {
  WriteToken(os, binary, elem);
}

template <std::size_t E, class T>
struct TupleWriterImpl {
  static constexpr size_t SIZE = std::tuple_size<T>::value;
  static void Write(std::ostream& os, bool binary, const T& tup) {
    constexpr std::size_t j = SIZE - E;
    /* WriteBasicType(os, binary, std::get<j>(tup)); */
    WriteTupleElement(os, binary, std::get<j>(tup));
    TupleWriterImpl<E - 1, T>::Write(os, binary, tup);
  }
};

template <class T>
struct TupleWriterImpl<0, T> {
  static void Write(std::ostream&, bool, const T&) { }
};

template <class T>
void WriteTuple(std::ostream &os, bool binary, const T& tup) {
  constexpr size_t SIZE = std::tuple_size<T>::value;
  TupleWriterImpl<SIZE, T>::Write(os, binary, tup);
}


// Utility classes to read an arbitrary tuple from an input stream

template <std::size_t E, class T>
struct TupleReaderImpl {
  static constexpr size_t SIZE = std::tuple_size<T>::value;
  static bool Read(std::istream& is, bool binary, T* tup) {
    constexpr std::size_t j = SIZE - E;
    const int i = is.peek();
    if (i == -1) {
      KALDI_WARN << "Unexpected EOF";
      return false;
    }
    // Some checks for a text input
    if (!binary) {
      if (static_cast<char>(i) == '\n') {
        KALDI_WARN << "Unexpected newline, reading vector< tuple<?> >; got "
                   << j << " elements, expected " << SIZE << ".";
        return false;
      } else if (static_cast<char>(i) == ';') {
        KALDI_WARN << "Wrong input format, reading vector< tuple<?> >; got "
                   << j << " elements, expected " << SIZE << ".";
        return false;
      }
    }
    ReadTuple(is, binary, &std::get<j>(tup));
    return TupleReaderImpl<E - 1, T>::Read(is, binary, tup);
  }
};

template <class T>
struct TupleReaderImpl<0, T> {
  static constexpr size_t SIZE = std::tuple_size<T>::value;
  static bool Read(std::istream& is, bool binary, T* /* unused */) {
    if (binary) {
      return true;
    } else {
      const int i = is.peek();
      if (i == -1) {
        KALDI_WARN << "Unexpected EOF";
        return false;
      } else if (static_cast<char>(i) == '\n') {
        return true;
      } else if (static_cast<char>(i) == ';') {
        return true;
      } else {
        KALDI_WARN << "Wrong input format, reading vector< tuple<?> >; got "
                   << "more than " << SIZE << " expected elements.";
        return false;
      }
    }
  }
};

template <class T>
bool ReadTuple(std::istream& is, bool binary, T* tup) {
  constexpr size_t SIZE = std::tuple_size<T>::value;
  return TupleReaderImpl<SIZE, T>::Read(is, binary, tup);
}

} // namespace internal

/// BasicTupleVectorHolder is a Holder for a vector of tuples of
/// a basic type, e.g. std::vector<std::tuple<int32, float> >.
/// Note: a basic type is defined as a type for which ReadBasicType
/// and WriteBasicType are implemented, i.e. integer and floating
/// types, and bool.
template<class ... BasicTypes> class BasicTupleVectorHolder {
 public:
  typedef std::tuple<BasicTypes...> Tuple;
  typedef std::vector< std::tuple<BasicTypes...> > T;
  static constexpr size_t tuple_size = std::tuple_size<Tuple>::value;

  BasicTupleVectorHolder() { }

  static bool Write(std::ostream &os, bool binary, const T &t) {
    InitKaldiOutputStream(os, binary);  // Puts binary header if binary mode.
    try {
      if (binary) {  // need to write the size, in binary mode.
        KALDI_ASSERT(static_cast<size_t>(static_cast<int32>(t.size())) ==
                     t.size());
        // Or this Write routine cannot handle such a large vector.
        // use int32 because it's fixed size regardless of compilation.
        // change to int64 (plus in Read function) if this becomes a problem.
        WriteBasicType(os, binary, static_cast<int32>(t.size()));
        for (typename T::const_iterator iter = t.begin();
             iter != t.end(); ++iter) {
          internal::WriteTuple(os, binary, *iter);
        }
      } else {  // text mode...
        // In text mode, we write out something like (for integers):
        // "1 2.2 ; 4 5.1 ; 6 -1.2 ; 8 3.2 \n"
        // where the semicolon is a separator, not a terminator.
        for (typename T::const_iterator iter = t.begin(); iter != t.end();) {
          internal::WriteTuple(os, binary, *iter);
          ++iter;
          if (iter != t.end())
            os << "; ";
        }
        os << '\n';
      }
      return os.good();
    } catch(const std::exception &e) {
      KALDI_WARN << "Exception caught writing Table object. " << e.what();
      return false;  // Write failure.
    }
    return false;
  }

  void Clear() { t_.clear(); }

  // Reads into the holder.
  bool Read(std::istream &is) {
    t_.clear();
    bool is_binary;
    if (!InitKaldiInputStream(is, &is_binary)) {
      KALDI_WARN << "Reading Table object [integer type], failed reading binary"
          " header\n";
      return false;
    }
    if (!is_binary) {
      // In text mode, we terminate with newline.
      try {  // catching errors from ReadBasicType..
        while (1) {
          const int i = is.peek();
          if (static_cast<char>(i) == '\n') {
            is.get();
            return true;
          } if (static_cast<char>(i) == ';') {
            is.get();
          } else if (static_cast<char>(i) == ' ') {
            is.get();
          } else {
            std::tuple<BasicTypes...> v;
            if (!internal::ReadTuple(is, is_binary, &v))
              return false;
            t_.push_back(v);
          }
        }
      } catch(const std::exception &e) {
        KALDI_WARN << "BasicTupleVectorHolder::Read, read error. " << e.what();
        return false;
      }
    } else {  // binary mode.
      size_t filepos = is.tellg();
      try {
        int32 size;
        ReadBasicType(is, true, &size);
        t_.resize(size);
        for (typename T::iterator iter = t_.begin(); iter != t_.end(); ++iter) {
          if (!internal::ReadTuple(is, is_binary, &(*iter)))
            return false;
        }
        return true;
      } catch(...) {
        KALDI_WARN << "BasicVectorHolder::Read, read error or unexpected data"
            " at archive entry beginning at file position " << filepos;
        return false;
      }
    }
    return false;
  }

  // Objects read/written with the Kaldi I/O functions always have the stream
  // open in binary mode for reading.
  static bool IsReadInBinary() { return true; }

  const T &Value() const {  return t_; }

  void Swap(BasicTupleVectorHolder<BasicTypes...> *other) {
    t_.swap(other->t_);
  }

  bool ExtractRange(const BasicTupleVectorHolder<BasicTypes...> &other,
                    const std::string &range) {
    KALDI_ERR << "ExtractRange is not defined for this type of holder.";
    return false;
  }

  ~BasicTupleVectorHolder() { }
 private:
  KALDI_DISALLOW_COPY_AND_ASSIGN(BasicTupleVectorHolder);

  T t_;
};


}  // namespace kaldi

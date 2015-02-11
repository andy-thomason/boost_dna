
#ifndef _SEED_H_
#define _SEED_H_

#include <ostream>
#include <cstdint>

//! This class contains a fixed-length DNA sequence using only the letters ACG and T.
//!
//! The intention of this class is to provide high performance low-level operators for
//! short DNA sequences. If longer, more complex sequences are required, refer to the
//! other classes.
//!
template<class _Base> class seed {
public:
  typedef _Base value_type;
  static const int num_bases = sizeof(value_type) * 4;

  seed(value_type value=0) : value_(value) {
  }

  template<class _FwdIter> seed(_FwdIter begin, _FwdIter end) {
    value_type value = 0;
    int seed_size = 0;
    for (_FwdIter p = begin; p != end && seed_size < num_bases; ++p) {
      value = value * 4 + (
        *p == 'C' ? 1 :
        *p == 'G' ? 2 :
        *p == 'T' ? 3 :
        0
      );
      seed_size++;
    }
    value_ = value << (num_bases - seed_size) * 2;
  }

  seed(const uint64_t *bp, uint64_t offset) {
    uint64_t sh = ((offset%32)*2);
    size_t idx = offset / 32;
    value_ = (bp[idx] << sh) | (bp[idx+1] >> (64-sh));
  }

  seed upper(int seed_size) const {
    return seed(value_ | ~(value_type)0 >> (seed_size*2));
  }

  //! The raw data value of this seed.
  value_type value() const {
    return value_;
  }

  //! 
  bool operator < (const seed &rhs) { return value_ < rhs.value_; }
  bool operator == (const seed &rhs) { return value_ == rhs.value_; }
  bool operator != (const seed &rhs) { return value_ != rhs.value_; }

  //! Code distance to another seed.
  int distance(const seed &rhs) const {
    value_type x = value_ ^ rhs.value_;
    x = (x & 0x5555555555555555) | ((x >> 1) & 0x5555555555555555);
    return __popcnt64(x);
  }

private:
  _Base value_;
};

typedef seed<std::uint32_t> seed32_t;
typedef seed<std::uint64_t> seed64_t;

template<class _Base> std::ostream &operator<<(std::ostream &os, const seed<_Base> &seed) {
  _Base value = seed.value();

  const int seed_size = seed.num_bases;
  char tmp[seed_size+1];
  for (int i = 0; i != seed_size; ++i) {
    tmp[i] = "ACGT"[value >> (seed_size*2-2)];
    value <<= 2;
  }
  tmp[seed_size] = 0;
  return os << tmp;
}

#endif


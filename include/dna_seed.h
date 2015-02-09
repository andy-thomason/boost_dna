
#ifndef DNA_SEED_INCLUDED
#define DNA_SEED_INCLUDED

#include <ostream>

//! This class contains a fixed-length DNA sequence using only the letters ACG and T.
//!
//! The intention of this class is to provide high performance low-level operators for
//! short DNA sequences. If longer, more complex sequences are required, refer to the
//! other classes.
//!
template<class _Base> class dna_seed {
  typedef _Base value_type;
public:
  dna_seed(value_type value) : value_(value) {
  }

  //! The raw data value of this seed.
  value_type value() const {
    return value_;
  }

  //! Code distance to another seed.
  int distance(const dna_seed &rhs) const {
    value_type xor = value ^ rhs.value;
    xor = (xor & 0x5555555555555555) | ((xor >> 1) & 0x5555555555555555);
    return __popcnt64(xor);
  }

private:
  _Base value_;
};

template<class _Base> std::ostream &operator<<(std::ostream &os, const dna_seed<_Base> &seed) {
  _Base value = seed.value;

  const int seed_size = sizeof(_Base) * 4;
  char tmp[seed_size+1];
  for (int i = 0; i != seed_size; ++i) {
    tmp[i] = "ACGT"[value >> (seed_size*2-2)];
    value <<= 2;
  }
  tmp[seed_size] = 0;
  return os;
}

#endif


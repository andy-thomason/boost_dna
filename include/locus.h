
#ifndef _LOCUS_H_
#define _LOCUS_H_

#include <ostream>

//! This class represents a location within 
//!
//!
template<class _Base> class locus {
  typedef _Base value_type;
public:
  locus(value_type value=0) : value_(value) {
  }

  //! The raw data value of this seed.
  value_type value() const {
    return value_;
  }

  bool operator<(const locus &rhs) const { return value_ < rhs.value_; }
  bool operator==(const locus &rhs) const { return value_ == rhs.value_; }
  bool operator!=(const locus &rhs) const { return value_ != rhs.value_; }

  int64_t operator-(const locus &rhs) const { return (int64_t)value_ - (int64_t)rhs.value_; }
  locus operator+(int64_t val) const { return locus(value_ + val); }

private:
  _Base value_;
};

typedef locus<uint64_t> locus64_t;
typedef locus<uint32_t> locus32_t;

template<class _Base> std::ostream &operator<<(std::ostream &os, const locus<_Base> &locus) {
  os << "[" << hex(locus.value()) << "]";
  return os;
}

#endif


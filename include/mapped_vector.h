
#ifndef _MAPPED_VECTOR_H_
#define _MAPPED_VECTOR_H_

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

template <class _Ty> class mapped_vector {
public:
  typedef _Ty value_type;

  mapped_vector() {
  }

  const value_type *data() const {
    return begin_;
  }

  size_t size() const {
    return end_ - begin_;
  }

  void map(size_t size, const char *begin) {
    begin_ = (const _Ty*)begin;
    end_ = begin_ + size;
  }

  const value_type &operator[](size_t index) const {
    return begin_[index];
  }
private:
  const _Ty *begin_;
  const _Ty *end_;
};

#endif


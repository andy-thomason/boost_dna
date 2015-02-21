
#ifndef _MAPPED_VECTOR_H_
#define _MAPPED_VECTOR_H_

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

template <class _Ty> class mapped_vector {
  boost::interprocess::file_mapping *mapping;
  boost::interprocess::mapped_region region;
  size_t size_;
  size_t capacity_;
  boost::interprocess::offset_t offset_;
public:
  typedef _Ty value_type;

  mapped_vector(
    boost::interprocess::file_mapping *mapping,
    boost::interprocess::offset_t offset,
    std::size_t size
  ) : mapping(mapping), size_(size), capacity_(size), offset_(offset) {
    region = boost::interprocess::mapped_region(
      *mapping, mapping->get_mode(), offset_, capacity_ * sizeof(value_type)
    );
  }

  mapped_vector() = delete;
  mapped_vector(const mapped_vector &) = delete;

  mapped_vector(const mapped_vector &&) {
  }

  value_type *data() { return (value_type *)region.get_address(); }
  size_t size() const { return size_; }
  size_t capacity() const { return capacity_; }
  const value_type *data() const { return (const value_type *)region.get_address(); }
  value_type &operator[](size_t index) { return data()[index]; }

  void push_back(const value_type &val) {
    size_t old_size = size_;
    resize(size_ + 1);
    data()[old_size] = val;
  }

  void resize(size_t new_size) {
    if (new_size > capacity_) {
      reserve(new_size + 0x100000 / sizeof(value_type));
    }
    size_ = new_size;
  }

  void reserve(size_t new_capacity) {
    if (new_capacity > capacity_) {
      capacity_ = new_capacity;
      region = boost::interprocess::mapped_region(
        mapping, mapping.get_mode(), offset_, capacity_ * sizeof(value_type)
      );
    }
  }
};

#endif


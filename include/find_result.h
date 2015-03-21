
#ifndef _FIND_RESULT_H_
#define _FIND_RESULT_H_

#include <ostream>
#include <string>
#include <algorithm>

class find_result {
  locus64_t locus_;
  std::string dna_;
  int errors_;
  std::string chromosome_;
  uint64_t offset_;
public:
  find_result() {
  }

  find_result(
    locus64_t locus,
    std::string dna,
    int errors,
    std::string chromosome,
    uint64_t offset
  ) : locus_(locus), dna_(dna), errors_(errors), chromosome_(chromosome), offset_(offset) {
  }

  locus64_t locus() const { return locus_; }
  std::string dna() const { return dna_; }
  int errors() const { return errors_; }
  std::string chromosome() const { return chromosome_; }
  int64_t offset() const { return offset_; }
};

std::ostream &operator <<(std::ostream &os, const find_result &r) {
  return os << r.chromosome() << " " << r.offset() << " " << r.errors() << " " << r.dna();
}

int dna_distance(const std::string &a, const std::string &b) {
  int errors = 0;
  size_t min_size = std::min(a.size(), b.size());
  for (size_t i = 0; i != min_size; ++i) {
    if (a[i] != b[i]) {
      errors += a[i] != 'N' && b[i] != 'N';
    }
  }
  return errors + std::abs((int)a.size() - (int)b.size());
}

#endif


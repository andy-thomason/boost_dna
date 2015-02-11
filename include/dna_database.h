
#ifndef _DNA_DATABASE_H_
#define _DNA_DATABASE_H_

#include "seed.h"
#include "locus.h"

#include <vector>
#include <cstdint>
#include <string>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cassert>

class dna_database {
  typedef std::vector<uint8_t> storage_t;
  typedef std::vector<uint8_t> aux_buffer_t;

  std::vector<uint64_t> bp_index;
  std::vector<uint64_t> bp;
  std::vector<uint64_t> aux_index;
  std::vector<uint8_t> aux;
  std::vector<uint8_t> D_data;
  uint64_t offset = 0;

  typedef std::pair<std::vector<uint8_t>::const_iterator, std::vector<uint8_t>::const_iterator> aux_result_t;

public:
  dna_database() {
  }

  void operator=(dna_database &rhs) = delete;
  dna_database(dna_database &rhs) = delete;


  class item {
  public:
    item(const dna_database *ptr, uint64_t index) : ptr(ptr), index(index) {
    }

    template <class _Fn> void for_each_forward_seed(_Fn fn, int seed_size) const {
      ptr->for_each_forward_seed(fn, index, seed_size);
    }

    void operator++() { index++; }
    void operator++(int) { index++; }
    bool operator!=(const item &rhs) { return index != rhs.index; }

    uint64_t size() {
      return ptr->get_size(index);
    }

    std::string &name(std::string &tmp) const {
      return ptr->get_name(tmp, index);
    }

    std::string &dna(std::string &tmp, bool newlines) const {
      return ptr->get_dna(tmp, index, newlines);
    }

    std::string &phred(std::string &tmp) const {
      return ptr->get_phred(tmp, index);
    }
  private:
    friend class dna_database;
    friend class iterator;
    const dna_database *ptr;
    uint64_t index;
  };

  class iterator : public item {
  public:
    iterator(const dna_database *ptr, uint64_t index) : item(ptr, index) {
    }

    iterator &operator++() { index++; return *this; }
    iterator &operator++(int) { index++; return *this; }
    bool operator!=(const iterator &rhs) { return index != rhs.index; }

    const item &operator *() const { return *this; }
    const item *operator->() const { return this; }
  };

  iterator begin() const {
    return iterator(this, 0);
  }

  iterator end() const {
    return iterator(this, bp_index.size() ? bp_index.size()-1 : 0);
  }

  std::ostream &dump(std::ostream &os) const {
    std::string tmp;
    if (true) {
      for (iterator p = begin(); p != end(); ++p) {
        os << ">" << p.name(tmp) << "\n";
        os << p.dna(tmp, true) << "\n";
      }
    } else {
      for (iterator p = begin(); p != end(); ++p) {
        os << "@" << p.name(tmp) << "\n";
        os << p.dna(tmp, false) << "\n";
        os << "+\n";
        //os << p.phred(tmp) << "\n";
      }
    }
    return os;
  }

  void begin_section() {
    if (!bp_index.size() || bp_index.back() != offset || aux_index.back() != aux.size()) {
      bp_index.push_back(offset);
      aux_index.push_back(aux.size());
    }
  }

  void add_name(const char *b, const char *e) {
    size_t len = e - b;
    aux.push_back((uint8_t)'N');
    push_vlq(aux, len);
    aux.insert(aux.end(), b, e);
  }

  void add_dna(const char *b, const char *e, const uint8_t *translate) {
    uint64_t bp_acc = 0;
    if (bp.size()*32 > offset) {
      bp_acc = bp.back() >> ((32-offset%32)*2);
      bp.pop_back();
    }

    uint8_t any_non_actg = 0;
    for (const char *q = b; q != e; ++q) {
      uint8_t t = translate[*q];
      if (!(t & 8)) {
        any_non_actg |= t;
        bp_acc = bp_acc * 4 + (t & 3);
        offset++;
        if ((offset & 31) == 0) {
          bp.push_back(bp_acc);
        }
      }
    }

    if (offset % 32) {
      bp.push_back(bp_acc << ((32-(offset % 32))*2));
    }

    if (any_non_actg & 4) {
      D_data.resize(0);
      for (const char *q = b; q != e;) {
        char c = *q++;
        uint8_t t = translate[c];
        size_t len = 1;
        if (!(t & 8)) {
          if (t & 4) {
            while (q != e) {
              char c2 = *q;
              uint8_t t = translate[c2];
              if (!(t & 8)) {
                if (c2 != c) break;
                len++;
              }
              ++q;
            }
            //std::cout << len << " x " << c << "\n";
            D_data.push_back(c);
            push_vlq(D_data, len);
          } else {
            while (q != e) {
              char c2 = *q;
              uint8_t t = translate[c2];
              if (!(t & 8)) {
                if (t & 4) break;
                len++;
              }
              ++q;
            }
            //std::cout << len << " x actg\n";
            D_data.push_back(0);
            push_vlq(D_data, len);
          }
        }
      }
      aux.push_back((uint8_t)'D');
      push_vlq(aux, D_data.size());
      aux.insert(aux.end(), D_data.begin(), D_data.end());
    }
  }

  uint64_t add_phred(const char *b, const char *e) {
    size_t len = e - b;
    aux.push_back((uint8_t)'P');
    push_vlq(aux, len);
    aux.insert(aux.end(), b, e);
  }

  void end_section() {
    aux.push_back(0);
    bp_index.push_back(offset);
    aux_index.push_back(aux.size());
  }

  const uint64_t *get_bp() const { return bp.data(); }

private:
  std::string &get_name(std::string &tmp, uint64_t index) const {
    aux_result_t ar = get_aux_data(index, 'N');
    return tmp.assign(ar.first, ar.second);
  }

  std::string &get_short_name(std::string &tmp, uint64_t index) const {
    aux_result_t ar = get_aux_data(index, 'N');
    auto p = ar.first;
    for (; p != ar.second && *p != ' '; ++p) {
    }
    return tmp.assign(ar.first, p);
  }

  std::string &get_dna(std::string &tmp, uint64_t index, bool newlines=false) const {
    locus64_t begin_offset = bp_index[index];
    locus64_t end_offset = bp_index[index+1];
    return get_dna_as_text(tmp, index, begin_offset, end_offset, newlines);
  }

  std::string &get_dna_as_text(std::string &tmp, uint64_t index, locus64_t begin_offset, locus64_t end_offset, bool newlines=false) const {
    tmp.resize((end_offset - begin_offset) + (newlines ? (end_offset - begin_offset)/60 : 0));

    size_t i = 0;
    uint64_t offset = bp_index[index];
    const uint64_t *p = &bp[offset/32];
    uint64_t p0 = p[0];
    uint64_t next_newline = newlines ? 60 : ~(uint64_t)0;

    aux_result_t ar = get_aux_data(index, 'D');
    if (ar.first != ar.second) {
      for (auto dp = ar.first; dp != ar.second && offset != end_offset.value(); ) {
        char c = *dp++;
        size_t len = 0;
        while (*dp & 0x80) {
          len |= (*dp++ & 0x7f);
          len <<= 7;
        }
        len |= (*dp++ & 0x7f);
        uint64_t run_end = std::min(offset + len, end_offset.value());
        if (c) {
          while (offset < run_end) {
            unsigned sh = (offset&31)*2;
            if (offset >= begin_offset.value()) {
              tmp[i++] = c;
              if (--next_newline == 0) { next_newline = 60; tmp[i++] = '\n'; }
            }
            if (++offset % 32 == 0) {
              p0 = *++p;
            }
          }
        } else {
          while (offset < run_end ) {
            unsigned sh = (offset&31)*2;
            if (offset >= begin_offset.value()) {
              uint64_t value = (p0 << sh >> 62) & 3;
              tmp[i++] = "ACGT"[value];
              if (--next_newline == 0) { next_newline = 60; tmp[i++] = '\n'; }
            }
            if (++offset % 32 == 0) {
              p0 = *++p;
            }
          }
        }
      }
    } else {
      while (offset != end_offset.value() ) {
        unsigned sh = (offset&31)*2;
        if (offset >= begin_offset.value()) {
          uint64_t value = (p0 << sh >> 62) & 3;
          tmp[i++] = "ACGT"[value];
          if (--next_newline == 0) { next_newline = 60; tmp[i++] = '\n'; }
        }
        if (++offset % 32 == 0) {
          p0 = *++p;
        }
      }
    }
    tmp.resize(i);
    return tmp;
  }

  std::string &get_phred(std::string &tmp, uint64_t index) const {
    aux_result_t ar = get_aux_data(index, 'P');
    tmp.assign(ar.first, ar.second);
    return tmp;
  }

  uint64_t get_size(uint64_t index) const {
    uint64_t begin_offset = bp_index[index];
    uint64_t end_offset = bp_index[index+1];
    return end_offset - begin_offset;
  }

  template <class _Fn> void for_each_forward_seed(_Fn fn, uint64_t index, int seed_size) const {
    uint64_t begin_offset = bp_index[index];
    uint64_t end_offset = bp_index[index+1] - seed_size + 1;

    uint64_t offset = begin_offset;
    const uint64_t *p = bp.data() + offset/32;
    uint64_t p0 = p[0];
    uint64_t p1 = p[1];
    int shift = (offset % 32) * 2;
    if (shift) { p0 = p0 << shift; p0 |= p1 >> (64-shift); p1 <<= shift; }
    uint64_t mask = ~(uint64_t)0 << ((32-seed_size)*2);

    while (offset < end_offset) {
      fn(offset, p0 & mask);
      p0 <<= 2; p0 |= p1 >> 62; p1 <<= 2;
      if (++offset % 32 == 0) {
        ++p;
        p0 = p[0];
        p1 = p[1];
      }
    }
  }

  // get (compressed) aux data for a section.
  aux_result_t get_aux_data(uint64_t index, uint16_t key) const {
    if (cur_aux != index) {
      cur_aux = index;
      const uint8_t *aux_data = aux.data();
      const uint8_t *begin = aux_data + aux_index[index];
      const uint8_t *end = aux_data + aux_index[index+1];
      // todo: decompress into aux_buffer
      aux_buffer.assign(begin, end);
    }

    auto e = aux_buffer.end();
    for (auto p = aux_buffer.begin(); *p && p != e; ) {
      bool matched = (p[0] & 0x80 ? ((p[0]&0x7f) * 256 + p[1]) : p[0]) == key;
      p += (p[0] & 0x80) ? 2 : 1;
      size_t len = 0;
      while (*p & 0x80) {
        len |= (*p++ & 0x7f);
        len <<= 7;
      }
      len |= *p++;
      if (matched) {
        return std::make_pair(p, p+len);
      }
      p += len;
    }
    return std::make_pair(e, e);
  }

  static void push_vlq(std::vector<uint8_t> &dest, size_t val) {
    if (val < 128) {
      dest.push_back((uint8_t)val);
    } else if (val < 128 * 128) {
      dest.push_back((uint8_t)(val >> 7) | 0x80);
      dest.push_back((uint8_t)(val) & 0x7f);
    } else if (val < 128 * 128 * 128) {
      dest.push_back((uint8_t)(val >> 14) | 0x80);
      dest.push_back((uint8_t)(val >> 7) | 0x80);
      dest.push_back((uint8_t)(val) & 0x7f);
    } else if (val < 128 * 128 * 128 * 128) {
      dest.push_back((uint8_t)(val >> 21) | 0x80);
      dest.push_back((uint8_t)(val >> 14) | 0x80);
      dest.push_back((uint8_t)(val >> 7) | 0x80);
      dest.push_back((uint8_t)(val) & 0x7f);
    }
  }

  size_t find_index(locus64_t locus) const {
    const uint64_t *b = bp_index.data();
    const uint64_t *e = b + bp_index.size();
    const uint64_t *p = std::upper_bound(b, e, locus.value());
    return p == b ? 0 : p - b - 1;
  }

  mutable aux_buffer_t aux_buffer;
  mutable uint64_t cur_aux = ~(uint64_t)0;

  friend class parser;
  friend class suffix_array;
};

std::ostream &operator<<(std::ostream &os, const dna_database &r) {
  return r.dump(os);
}

class parser {
  dna_database *dna;
  uint8_t translate[256];

  void err(const char *msg) {
  }

public:
  parser() {
    memset(translate, 4, sizeof(translate));
    translate['\r'] = 8;
    translate['\n'] = 8;
    translate['A'] = 0;
    translate['C'] = 1;
    translate['G'] = 2;
    translate['T'] = 3;
  }

  const char *add_fasta(dna_database *dna, const char *begin, const char *end, size_t max_chromosomes=~(size_t)0) {
    this->dna = dna;

    std::vector<uint8_t> &aux = dna->aux;

    const char * p = begin;
    for (size_t num_chromosomes = 0; p != end && *p == '>' && num_chromosomes < max_chromosomes; ++num_chromosomes) {
      dna->begin_section();

      const char *b = p;
      for (; p != end && *p != '\n'; ++p) {
      }

      dna->add_name(b+1, p);

      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      b = p;
      while (p != end && *p != '>') ++p;
      dna->add_dna(b, p, translate);

      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      dna->end_section();
    }

    return p;
  }

  void add_fastq(dna_database *dna, const char * begin, const char * end) {
    this->dna = dna;

    std::vector<uint8_t> &aux = dna->aux;

    const char *p = begin;
    while (p != end) {
      dna->begin_section();

      if (*p != '@') {
        err("expected @ in fq file");
      }

      const char *b = p;
      for (; p != end && *p != '\n'; ++p) {
      }
      dna->add_name(b+1, p);
      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      b = p;
      for (; p != end && *p != '\n'; ++p) {
      }
      dna->add_dna(b, p, translate);
      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      b = p;
      for (; p != end && *p != '\n'; ++p) {
      }
      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      b = p;
      for (; p != end && *p != '\n'; ++p) {
      }

      dna->add_phred(b, p);

      while (p != end && (*p == '\n' || *p == '\r')) ++p;
      dna->end_section();
    }
  }
};

#endif


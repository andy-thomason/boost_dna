
#ifndef _DNA_DATABASE_H_
#define _DNA_DATABASE_H_

#include "seed.h"
#include "locus.h"
#include "find_result.h"

#include <vector>
#include <cstdint>
#include <string>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cassert>

template <class bp_index_vector_type, class bp_vector_type, class aux_index_vector_type, class aux_vector_type> class dna_database {
  bp_index_vector_type bp_index;
  bp_vector_type bp;
  aux_index_vector_type aux_index;
  aux_vector_type aux;
  typedef size_t item_index_t;

  std::vector<uint8_t> D_data;
public:
  typedef typename bp_index_vector_type::value_type bp_index_type;
  typedef typename bp_vector_type::value_type bp_type;
  typedef typename aux_index_vector_type::value_type aux_index_type;
  typedef typename aux_vector_type::value_type aux_type;
  typedef std::pair<const aux_type*, const aux_type*> aux_result_t;

  dna_database() {
  }

  void operator=(dna_database &rhs) = delete;
  dna_database(dna_database &rhs) = delete;


  class item {
  public:
    item(const dna_database *ptr, item_index_t index) : ptr(ptr), index(index) {
    }

    template <class _Fn> void for_each_seed(_Fn fn, int seed_size) const {
      ptr->for_each_seed(fn, index, seed_size);
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
    item_index_t index;
  };

  class iterator : public item {
  public:
    iterator(const dna_database *ptr, item_index_t index) : item(ptr, index) {
    }

    iterator &operator++() { index++; return *this; }
    iterator &operator++(int) { ++index; return *this; }
    bool operator!=(const iterator &rhs) { return index != rhs.index; }

    const item &operator *() const { return *this; }
    const item *operator->() const { return this; }
  };

  iterator begin() const {
    return iterator(this, 0);
  }

  iterator end() const {
    return iterator(this, size());
  }

  template <class _Write> bool write(_Write &w) const {
    static const char sig[32] = "dna_database v1.0\n\x1a\x04";
    w.write(sig, sizeof(sig));

    size_t bp_index_size = bp_index.size();
    size_t bp_size = bp.size();
    size_t aux_index_size = aux_index.size();
    size_t aux_size = aux.size();
    w.write((const char*)&bp_index_size, sizeof(size_t));
    w.write((const char*)&bp_size, sizeof(size_t));
    w.write((const char*)&aux_index_size, sizeof(size_t));
    w.write((const char*)&aux_size, sizeof(size_t));

    w.write((const char*)bp_index.data(), bp_index_size * sizeof(bp_index_type));
    w.write((const char*)bp.data(), bp_size * sizeof(bp_type));
    w.write((const char*)aux_index.data(), aux_index_size * sizeof(aux_index_type));
    w.write((const char*)aux.data(), aux_size * sizeof(aux_type));
    return !w.fail();
  }

  template <class _Read> bool read(_Read &r) {
    char sig[32];
    r.read(sig, sizeof(sig));
    if (strcmp(sig, "dna_database v1.0\n\x1a\x04")) {
      return false;
    }

    size_t bp_index_size = 0;
    size_t bp_size = 0;
    size_t aux_index_size = 0;
    size_t aux_size = 0;

    r.read((char*)&bp_index_size, sizeof(size_t));
    r.read((char*)&bp_size, sizeof(size_t));
    r.read((char*)&aux_index_size, sizeof(size_t));
    r.read((char*)&aux_size, sizeof(size_t));

    bp_index.resize(bp_index_size);
    bp.resize(bp_size);
    aux_index.resize(aux_index_size);
    aux.resize(aux_size);

    r.read((char*)bp_index.data(), bp_index_size * sizeof(bp_index_type));
    r.read((char*)bp.data(), bp_size * sizeof(bp_type));
    r.read((char*)aux_index.data(), aux_index_size * sizeof(aux_index_type));
    r.read((char*)aux.data(), aux_size * sizeof(aux_type));
    return !r.fail();
  }

  void reserve(size_t bp_index_size, size_t bp_size, size_t aux_index_size, size_t aux_size) {
    bp_index.reserve(bp_index_size);
    bp.reserve(bp_size);
    aux_index.reserve(aux_index_size);
    bp.reserve(aux_size);
  }

  void begin_section() {
    if (!bp_index.size()) {
      bp_index.push_back((bp_index_type)0);
      aux_index.push_back((aux_index_type)0);
    }
    bp_index.push_back((bp_index_type)bp_index.back());
    aux_index.push_back((aux_index_type)aux_index.back());
  }

  void add_name(const char *b, const char *e) {
    size_t len = e - b;
    aux.push_back((uint8_t)'N');
    push_vlq(aux, len);
    aux.insert(aux.end(), b, e);
  }

  void add_dna(const char *b, const char *e) {
    // 0-3=ACGT 4=N etc. 8=ignore
    static const uint8_t translate[] = {
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8,
      8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    };
    uint64_t bp_acc = 0;
    bp_index_type offset = bp_index.back();
    if (bp.size()*32 > offset) {
      bp_acc = bp.back() >> ((32-offset%32)*2);
      bp.pop_back();
    }

    size_t extra = ((e - b) + (sizeof(bp_type) * 4) - 1 ) / (sizeof(bp_type) * 4);
    bp.reserve(bp.size() + extra);

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
    bp_index.back() = offset;

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
      aux_index.back() = (aux_index_type)aux.size();
    }
  }

  void add_phred(const char *b, const char *e) {
    size_t len = e - b;
    aux.push_back((uint8_t)'P');
    push_vlq(aux, len);
    aux.insert(aux.end(), b, e);
    aux_index.back() = (aux_index_type)aux.size();
  }

  void end_section() {
    aux.push_back(0);
    aux_index.back() = (aux_index_type)aux.size();
  }

  const bp_type *bp_data() const { return bp.data(); }
  const bp_index_type *bp_index_data() const { return bp_index.data(); }

  size_t size() const {
    return bp_index.size() ? bp_index.size()-1 : 0;
  }

  std::string &get_name(std::string &tmp, item_index_t index) const {
    aux_result_t ar = get_aux_data(index, 'N');
    return tmp.assign(ar.first, ar.second);
  }

  std::string &get_short_name(std::string &tmp, item_index_t index) const {
    aux_result_t ar = get_aux_data(index, 'N');
    auto p = ar.first;
    for (; p != ar.second && *p != ' '; ++p) {
    }
    return tmp.assign(ar.first, p);
  }

  std::string &get_dna(std::string &tmp, item_index_t index, bool newlines=false) const {
    locus64_t begin_offset = bp_index[index];
    locus64_t end_offset = bp_index[index+1];
    return get_dna_as_text(tmp, index, begin_offset, end_offset, newlines);
  }

  locus64_t get_dna_start(item_index_t index) const {
    return bp_index[index];
  }

  locus64_t get_dna_end(item_index_t index) const {
    return bp_index[index+1];
  }

  std::string &get_dna_as_text(std::string &tmp, locus64_t begin_offset, locus64_t end_offset, bool newlines=false) const {
    tmp.resize((end_offset - begin_offset) + (newlines ? (end_offset - begin_offset)/60 : 0));

    size_t i = 0;
    uint64_t next_newline = newlines ? 60 : ~(uint64_t)0;

    size_t index = find_index(begin_offset);
    aux_result_t ar = get_aux_data(index, 'D');
    if (0 && ar.first != ar.second) {
      uint64_t offset = bp_index[index];
      const bp_type *p = &bp[offset/32];
      uint64_t p0 = p[0];
      for (auto dp = ar.first; dp != ar.second && offset < end_offset.value(); ) {
        char c = *dp++;
        size_t len = 0;
        while (*dp & 0x80) {
          len |= (*dp++ & 0x7f);
          len <<= 7;
        }
        len |= (*dp++ & 0x7f);
        uint64_t run_end = std::min(offset + len, end_offset.value());
        if (run_end <= begin_offset.value()) {
          offset = run_end;
        } else {
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
            while (offset < run_end) {
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
      }
    } else {
      uint64_t offset = begin_offset.value();
      const bp_type *p = &bp[offset/32];
      uint64_t p0 = p[0];
      while (offset != end_offset.value() ) {
        unsigned sh = (offset&31)*2;
        uint64_t value = (p0 << sh >> 62) & 3;
        tmp[i++] = "ACGT"[value];
        if (--next_newline == 0) { next_newline = 60; tmp[i++] = '\n'; }
        if (++offset % 32 == 0) {
          p0 = *++p;
        }
      }
    }
    tmp.resize(i);
    return tmp;
  }

  std::string &get_phred(std::string &tmp, item_index_t index) const {
    aux_result_t ar = get_aux_data(index, 'P');
    tmp.assign(ar.first, ar.second);
    return tmp;
  }

  uint64_t get_size(item_index_t index) const {
    uint64_t begin_offset = bp_index[index];
    uint64_t end_offset = bp_index[index+1];
    return end_offset - begin_offset;
  }

  template <class _Fn> void for_each_seed(_Fn fn, item_index_t index, int seed_size) const {
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

  template <class _Dest, class _Filter> _Dest *copy_if(_Dest *__restrict dest, _Dest *dest_max, _Filter filter) __restrict {
    const bp_type *__restrict bp = bp_data();
    const bp_index_type *__restrict bp_index = bp_index_data();
    size_t max_idx = size();
    for (size_t idx = 0; idx != max_idx; ++idx) {
      bp_index_type begin_offset = bp_index[idx];
      bp_index_type end_offset = bp_index[idx+1] - 32 + 1;

      bp_index_type offset = begin_offset;
      const uint64_t *__restrict p = bp + offset/32;
      uint64_t p0 = p[0];
      uint64_t p1 = p[1];
      int shift = (offset % 32) * 2;
      if (shift) { p0 = p0 << shift; p0 |= p1 >> (64-shift); p1 <<= shift; }

      while (offset < end_offset) {
        bp_index_type next_offset = std::min((offset + 32) & -32, end_offset);
        while (offset < next_offset) {
          if (filter(seed64_t(p0))) {
            if (dest == dest_max) return dest;
            *dest++ = _Dest(p0, offset);
          }
          p0 <<= 2; p0 |= p1 >> 62; p1 <<= 2;
          ++offset;
        }
        ++p;
        p0 = p[0];
        p1 = p[1];
      }
    }
    return dest;
  }

  size_t find_index(locus64_t locus) const {
    const bp_index_type *b = bp_index.data();
    const bp_index_type *e = bp_index.data() + bp_index.size();
    const bp_index_type *p = std::upper_bound(b, e, locus.value());
    return p == b ? 0 : p - b - 1;
  }

  // brute force find
  bool find(std::vector<find_result> &result, const std::string &sequence, int max_distance, int max_results) {
    typedef std::pair<uint64_t, uint64_t> match_t;
    std::vector<match_t> matches(max_results ? max_results : 0x10000);
    seed64_t test_seed(sequence.begin(), sequence.end());
    int seed_size = std::min((int)sequence.size(), 32);
    uint64_t mask = ~(uint64_t)0 << (64-seed_size*2);
    match_t *last = copy_if(
      matches.data(), matches.data() + matches.size(),
      [=](seed64_t seed) {
        int distance = test_seed.distance(seed, mask);
        return distance <= max_distance;
      }
    );

    std::string text;
    result.resize(0);
    for (match_t *p = matches.data(); p != last; ++p) {
      //printf("%llx %llx\n", p->first, p->second);
      size_t index = find_index(p->second);
      locus64_t start = get_dna_start(index);
      std::string short_name;
      get_short_name(short_name, index);
      result.push_back(find_result{p->second, text, 0, short_name, p->second - start.value() + 1});
    }

    return true;
  }

private:
  // get (compressed) aux data for a section.
  aux_result_t get_aux_data(item_index_t index, uint16_t key) const {
    const aux_type *aux_buffer_begin = aux.data() + aux_index[index];
    const aux_type *aux_buffer_end = aux.data() + aux_index[index+1];

    for (const aux_type *p = aux_buffer_begin; *p && p != aux_buffer_end; ) {
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
    return std::make_pair(aux_buffer_end, aux_buffer_end);
  }

  template <class vector_type> static void push_vlq(vector_type &dest, size_t val) {
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

  mutable std::vector<uint8_t> aux_buffer;
  mutable uint64_t cur_aux = ~(uint64_t)0;
};

template <class... args> std::ostream &operator<<(std::ostream &os, const dna_database<args...> &r) {
  //return r.dump(os);
}

class parser {
  void err(const char *msg) {
  }

public:
  parser() {
  }

  template<class _Dna> const char *add_fasta(_Dna *dna, const char *begin, const char *end, size_t max_chromosomes=~(size_t)0) {
    dna->reserve(256, (end - begin) / 32, 256, 0x10000);

    const char * p = begin;
    for (size_t num_chromosomes = 0; p != end && *p == '>' && num_chromosomes < max_chromosomes; ++num_chromosomes) {
      dna->begin_section();

      const char *b = p;
      for (; p != end && *p != '\n'; ++p) {
      }

      dna->add_name(b+1, p);
      std::cout << std::string(b+1, p) << std::endl;

      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      b = p;
      while (p != end && *p != '>') ++p;
      dna->add_dna(b, p);

      while (p != end && (*p == '\n' || *p == '\r')) ++p;

      dna->end_section();
    }

    return p;
  }

  template<class _Dna> void add_fastq(_Dna *dna, const char * begin, const char * end) {
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


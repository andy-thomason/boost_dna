
#ifndef _SUFFIX_ARRAY_H_
#define _SUFFIX_ARRAY_H_

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
#include <thread>
#include <atomic>
#include <mutex>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <future>
#include <type_traits>

#include "mapped_vector.h"

template <class _Fn> static void par_for(int limit, _Fn f) {
  std::vector<std::future<bool>> threads(std::thread::hardware_concurrency());
  //std::vector<std::future<bool>> threads(1);

  std::atomic<int> sequence;
  for (auto &t : threads) {
    t = std::async([&](){
      int seq;
      while ((seq = sequence++) < limit) {
        f(seq);
      }
      return true;
    });
  }

  for (auto &t : threads) {
    t.get();
  }
}

template <class index_vector_type, class address_vector_type> class suffix_array {
  index_vector_type index;
  address_vector_type address;
  typedef typename index_vector_type::value_type index_type;
  typedef typename address_vector_type::value_type address_type;
public:
  suffix_array() {
  }

  template <class _Dna> suffix_array(_Dna &dna, int lg_index_size) {
    index.resize(((size_t)1 << lg_index_size)+1);

    enum { lg_num_seqs = 5 };
    enum { num_seqs = 1 << lg_num_seqs };
    size_t counts[num_seqs];
    std::fill(std::begin(counts), std::end(counts), 0);

    std::cout << "counting\n";
    size_t count = 0;
    int num_values = 0;
    std::string tmp;
    for (auto p = dna.begin(); p != dna.end(); ++p) {
      std::cout << p.name(tmp) << "\n";
      const uint64_t *bp = dna.bp_data();
      p->for_each_seed(
        [&](uint64_t offset, uint64_t seed) {
          //if (!seed) return;
          //std::cout << hex(seed) << " " << hex(seed64_t(bp, offset).value()) << "\n";
          if (seed) counts[seed >> (64-lg_num_seqs)]++;
        }, 32
      );
    }

    const uint64_t *bp = dna.bp_data();

    address.resize(std::accumulate(std::begin(counts), std::end(counts), (size_t)0));
    size_t starts[num_seqs];
    std::partial_sum(std::begin(counts), std::end(counts), starts);

    std::ofstream dbg("dbg.txt");

    par_for(num_seqs, [&](int seq) {
      std::cout << "seq=" << hex(seq) << " counts[seq]=" << counts[seq] << "\n";
      struct sort_t {
        uint64_t val;
        sort_t() {
        }
        sort_t(uint64_t seed, uint64_t offset) {
          _mm_stream_si64((int64_t*)&this->val, (seed >> (32-lg_num_seqs) << 32) | (offset & 0xffffffff));
        }
        bool operator<(sort_t rhs) { return val < rhs.val; }
      };
      //typedef std::pair<uint32_t, _Dna::bp_index_type> sort_t;
      long long t1 = __rdtsc();
      sort_t *sorter = new sort_t[counts[seq]];
      size_t sorter_size = counts[seq];
      long long t2 = __rdtsc();
      sort_t *last = dna.copy_if(
        sorter, sorter + sorter_size,
        [=](seed64_t seed) {
          if (seed == 0) return false;
          return (seed.value() >> (64 - lg_num_seqs)) == seq;
        }
      );
      printf("%lld %lld\n", last - sorter, sorter_size);
      assert(last - sorter == sorter_size);

      long long t3 = __rdtsc();
      std::sort(sorter, sorter + sorter_size);

      long long t4 = __rdtsc();
      size_t addr_idx = seq == 0 ? (size_t)0 : starts[seq-1];
      size_t first_idx_idx = (size_t)seq << (lg_index_size - lg_num_seqs);
      size_t last_idx_idx = (size_t)(seq+1) << (lg_index_size - lg_num_seqs);
      size_t idx_idx = first_idx_idx;
      uint64_t next_val = 0;
      int sh = 64 - lg_index_size + lg_num_seqs;

      for (size_t i = 0; i != sorter_size; ++i) {
        sort_t &p = sorter[i];
        while ((p.val >> sh) >= next_val) {
          index[idx_idx++] = (uint32_t)addr_idx;
          assert(idx_idx < last_idx_idx);
          next_val++;
        }
        //dbg << hex((uint8_t)(seq)) <<  " " << hex((uint32_t)(p.val >> sh)) << " " << hex((uint64_t)p.val) << " " << hex(seed64_t(bp, (uint32_t)p.val).value()) << "\n";
        address[addr_idx++] = (uint32_t)p.val;
      }
      assert(addr_idx == starts[seq]);
      //std::cout << hex(first_idx_idx) << "/" << hex(idx_idx) << "/" << hex(last_idx_idx) << "\n";
      while (idx_idx < last_idx_idx) {
        index[idx_idx++] = (uint32_t)addr_idx;
      }
      long long t5 = __rdtsc();
      printf("%lld %lld %lld %lld\n", t2-t1, t3-t2, t4-t3, t5-t4);
      delete [] sorter;
    });

    index.back() = (index_type)address.size();

    /*for (uint16_t idx_idx = 0; idx_idx != (1<<lg_index_size); ++idx_idx) {
      std::cout << hex(idx_idx) << " " << index[idx_idx] << std::endl;
      for (uint16_t addr_idx = index[idx_idx]; addr_idx != index[idx_idx+1]; ++addr_idx) {
        std::cout << address[addr_idx] << "/" << seed64_t(bp, address[addr_idx].value()) << "\n";
      }
      std::cout << "\n";
    }*/
  }

  template <class _Dna> bool verify(_Dna &dna) const {
    if (index.size() == 0) return false;

    const uint64_t *bp = dna.bp_data();
    size_t imax = index.size() - 1;
    for (size_t i = 0; i != imax; ++i) {
      index_type begin = index[i];
      index_type end = index[i+1];
      for (index_type j = begin; j != end; ++j) {
        address_type addr = address[j];
        seed64_t seed(bp, addr);
        std::cout << seed << "\n";
      }
    }
    return false;
  }

  void operator=(suffix_array &rhs) = delete;

  void operator=(suffix_array &&rhs) {
    index = std::move(rhs.index);
    address = std::move(rhs.address);
  }

  class iterator {
    size_t addr_idx;
    const suffix_array &suffix;
  public:
    iterator(const suffix_array &suffix, size_t addr_idx) : suffix(suffix), addr_idx(addr_idx) {
    }
  };

  iterator begin() const {
    return iterator(*this, 0);
  }

  iterator end() const {
    return iterator(*this, address.size());
  }

  template <class _Dna> bool find_seeds(_Dna &dna, std::vector<locus64_t> &result, seed64_t lower, seed64_t upper, int offset, int max_distance, int max_results) {
    int lg_index_size = (int)std::log2((int)index.size());
    size_t lower_idx = lower.value() >> (64-lg_index_size);
    size_t upper_idx = upper.value() >> (64-lg_index_size);

    const address_type *begin = address.data() + index[lower_idx];
    const address_type *end = address.data() + index[upper_idx+1];
    const uint64_t *bp = dna.bp_data();

    //std::cout << "find_seed " << lower << ".." << upper << " idx=" << lower_idx << "\n";

    const address_type *p = std::lower_bound(
      begin, end, lower, [=](address_type a, seed64_t b) {
        seed64_t seed_at_a(bp, a);
        //std::cout << seed_at_a << " < " << b << "\n";
        return seed_at_a < b;
      }
    );

    while (p < end) {
      seed64_t seed_at_p(bp, *p);
      if (!(seed_at_p < upper)) break;
      result.push_back(*p - offset);
      ++p;
    }

    return true;
  }

  template <class _Dna> bool find(_Dna &dna, std::vector<find_result> &result, const std::string &sequence, int max_distance, int max_results) {
    result.resize(0);
    int seq_size = (int)sequence.size();

    long long t1 = __rdtsc();

    // AAAAXAAAAA (10-1)/(1+1) = 4
    // AAXAAAXAAA (10-2)/(2+1) = 2
    // AXAAXAAXAA (10-3)/(3+1) = 1
    // AAAAAXAAAAA (11-1)/(1+1) = 5
    // AAAXAAAXAAA (11-2)/(2+1) = 3
    // AAXAAXAAXAA (11-3)/(3+1) = 2

    int seed_size = (seq_size - max_distance) / (max_distance+1);
    if (seed_size <= 0) return false;
    seed_size = std::min(seed_size, 32);
    const char *begin = sequence.data();
    std::vector<locus64_t> loci;
    for (int i = 0; i <= seq_size - seed_size; ++i) {
      seed64_t seed(begin + i, begin + i + seed_size);
      if (!find_seeds(dna, loci, seed, seed.upper(seed_size), i, 0, max_results)) {
        return false;
      }
    }

    long long t2 = __rdtsc();

    std::sort(loci.begin(), loci.end());

    long long t3 = __rdtsc();

    std::string text;
    for (auto p = loci.begin(); p != loci.end(); ) {
      auto q = p + 1;
      uint64_t index = dna.find_index(*p);
      dna.get_dna_as_text(text, index, *p, *p + seq_size, false);
      //std::cout << *p << " " << text << "\n";
      int errors = dna_distance(sequence, text);
      if (errors <= max_distance) {
        size_t index = dna.find_index(*p);
        locus64_t start = dna.get_dna_start(index);
        std::string short_name;
        dna.get_short_name(short_name, index);
        result.push_back(find_result{*p, text, errors, short_name, (uint64_t)(*p - start + 1)});
        if (max_results && loci.size() > max_results) {
          return false;
        }
      }
      for (; q != loci.end() && *q == *p; ++q) {
      }
      p = q;
    }

    long long t4 = __rdtsc();

    std::cout << t2-t1 << " " << t3-t2 << " " << t4-t3 << "\n";
    return true;
  }

  std::ostream &dump(std::ostream &os) {
    std::string tmp;
    return os;
  }

  template <class _Writeable> bool write(_Writeable &w) const {
    static const char sig[32] = "suffix_array v1.0\n\x1a\x04";
    w.write(sig, sizeof(sig));

    size_t index_size = index.size();
    size_t address_size = address.size();
    w.write((const char*)&index_size, sizeof(size_t));
    w.write((const char*)&address_size, sizeof(size_t));

    w.write((const char*)index.data(), index.size() * sizeof(index_type));
    w.write((const char*)address.data(), address.size() * sizeof(address_type));
    return !w.fail();
  }

  template <class _Readable> bool read(_Readable &r) {
    char sig[32];
    r.read(sig, sizeof(sig));
    if (strcmp(sig, "suffix_array v1.0\n\x1a\x04")) {
      return false;
    }

    size_t index_size = 0;
    size_t address_size = 0;

    r.read((char*)&index_size, sizeof(size_t));
    r.read((char*)&address_size, sizeof(size_t));

    index.resize(index_size);
    address.resize(address_size);

    char *index_data = (char*)index.data();
    r.read(index_data, index_size * sizeof(index_type));
    r.read((char*)address.data(), address.size() * sizeof(address_type));
    return !r.fail();
  }
};

#endif


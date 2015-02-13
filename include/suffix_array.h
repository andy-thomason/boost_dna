
#ifndef _SUFFIX_ARRAY_H_
#define _SUFFIX_ARRAY_H_

#include "dna_database.h"
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
#include <thread>
#include <atomic>
#include <mutex>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <future>

class suffix_array {
  dna_database &dna;
  int lg_index_size;

  std::vector<uint32_t> index;
  //std::vector<uint8_t> remainder;
  std::vector<locus32_t> addresses;

  template <class _Fn> static void par_for(int limit, _Fn f) {
    //for (int i = 0; i != limit; ++i) {
      //f(i);
    //}
    std::vector<std::future<bool>> threads(std::thread::hardware_concurrency());

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

  /*template <class... args_t> void log(std::string &str, args_t... args) {
    str.resize(sprintf(nullptr, args...));
    sprintf(str.data(), args...);
    std::cout << str;
  }*/

public:
  suffix_array(dna_database &dna, int lg_index_size) : dna(dna), lg_index_size(lg_index_size) {
    index.resize(((size_t)1 << lg_index_size)+1);

    const int lg_num_seqs = 4;
    const int num_seqs = 1 << lg_num_seqs;
    size_t counts[num_seqs];
    std::fill(std::begin(counts), std::end(counts), 0);

    par_for(num_seqs, [&](int seq) {
      size_t count = 0;
      int num_values = 0;
      for (auto p = dna.begin(); p != dna.end(); ++p) {
        int seed_size = 32;
        p->for_each_forward_seed(
          [&](uint64_t offset, uint64_t seed) {
            if ((seed >> (64-lg_num_seqs)) == seq) {
              count++;
            }
          }, seed_size
        );
      }
      counts[seq] = count;
    });

    const uint64_t *bp = dna.bp.data();

    addresses.resize(std::accumulate(std::begin(counts), std::end(counts), (size_t)0));
    size_t starts[num_seqs];
    std::partial_sum(std::begin(counts), std::end(counts), starts);

    par_for(num_seqs, [&](int seq) {
      std::vector<std::pair<uint64_t, uint64_t>> sorter;
      sorter.resize(counts[seq]);
      size_t count = 0;
      for (auto p = dna.begin(); p != dna.end(); ++p) {
        int seed_size = 32;
        p->for_each_forward_seed(
          [&](uint64_t offset, uint64_t seed) {
            if ((seed >> (64-lg_num_seqs)) == seq) {
              sorter[count++] = std::make_pair(seed, offset);
            }
          }, seed_size
        );
      }

      std::sort(sorter.begin(), sorter.end());

      size_t addr_idx = seq == 0 ? (size_t)0 : starts[seq-1];
      size_t idx_idx = (size_t)seq << (lg_index_size - lg_num_seqs);
      size_t last_idx_idx = (size_t)(seq+1) << (lg_index_size - lg_num_seqs);
      for (auto p : sorter) {
        while ((p.first >> (64-lg_index_size)) >= idx_idx) {
          index[idx_idx++] = (uint32_t)addr_idx;
        }
        //std::cout << hex(idx_idx) << " " << hex(p.first) << " " << hex(p.second) << " " << seed64_t(bp, p.second) << "\n";
        addresses[addr_idx++] = (uint32_t)p.second;
      }
      assert(addr_idx == starts[seq]);
      while (idx_idx < last_idx_idx) {
        index[idx_idx++] = (uint32_t)addr_idx;
      }
    });

    index.back() = (uint32_t)addresses.size();

    /*for (uint16_t idx_idx = 0; idx_idx != (1<<lg_index_size); ++idx_idx) {
      std::cout << hex(idx_idx) << " " << index[idx_idx] << std::endl;
      for (uint16_t addr_idx = index[idx_idx]; addr_idx != index[idx_idx+1]; ++addr_idx) {
        std::cout << addresses[addr_idx] << "/" << seed64_t(bp, addresses[addr_idx].value()) << "\n";
      }
      std::cout << "\n";
    }*/
  }

  void operator=(suffix_array &rhs) = delete;
  suffix_array(dna_database &rhs) = delete;

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
    return iterator(*this, addresses.size());
  }

  bool find_seeds(std::vector<locus64_t> &result, seed64_t lower, seed64_t upper, int offset, int max_distance, int max_results) {
    size_t lower_idx = lower.value() >> (64-lg_index_size);
    size_t upper_idx = upper.value() >> (64-lg_index_size);

    const locus32_t *begin = addresses.data() + index[lower_idx];
    const locus32_t *end = addresses.data() + index[upper_idx+1];
    const uint64_t *bp = dna.bp.data();

    //std::cout << "find_seed " << lower << ".." << upper << " idx=" << lower_idx << "\n";

    const locus32_t *p = std::lower_bound(
      begin, end, lower, [=](locus32_t a, seed64_t b) {
        seed64_t seed_at_a(bp, a.value());
        return seed_at_a < b;
      }
    );

    while (p < end) {
      seed64_t seed_at_p(bp, p->value());
      if (!(seed_at_p < upper)) break;
      result.push_back(p->value() - offset);
      ++p;
    }

    return true;
  }

  struct find_result {
    locus64_t locus_;
    std::string dna;
    int errors;
    std::string chromosome;
    int64_t offset;
  };

  bool find(std::vector<find_result> &result, const std::string &sequence, int max_distance, int max_results) {
    result.resize(0);
    int seq_size = (int)sequence.size();

    // AAAAXAAAAA (10-1)/(1+1) = 4
    // AAXAAAXAAA (10-2)/(2+1) = 2
    // AXAAXAAXAA (10-3)/(3+1) = 1
    // AAAAAXAAAAA (11-1)/(1+1) = 5
    // AAAXAAAXAAA (11-2)/(2+1) = 3
    // AAXAAXAAXAA (11-3)/(3+1) = 2

    int seed_size = max_distance == 0 ? seq_size : (seq_size - max_distance) / max_distance;
    if (seed_size <= 0) return false;
    seed_size = std::min(seed_size, 32);
    const char *begin = sequence.data();
    std::vector<locus64_t> loci;
    for (int i = 0; i <= seq_size - seed_size; ++i) {
      seed64_t seed(begin + i, begin + i + seed_size);
      if (!find_seeds(loci, seed, seed.upper(seed_size), i, 0, max_results)) {
        return false;
      }
    }

    std::sort(loci.begin(), loci.end());

    std::string text;
    for (auto p = loci.begin(); p != loci.end(); ) {
      auto q = p + 1;
      uint64_t index = dna.find_index(*p);
      dna.get_dna_as_text(text, index, *p, *p + seq_size, false);
      //std::cout << *p << " " << text << "\n";
      int errors = 0;
      for (size_t i = 0; i != seq_size; ++i) {
        if (sequence[i] != text[i] || sequence[i] == 'N' || text[i] == 'N') {
          errors++;
          if (errors > max_distance) break;
        }
      }
      if (errors <= max_distance) {
        size_t index = dna.find_index(*p);
        locus64_t start = dna.get_dna_start(index);
        std::string short_name;
        dna.get_short_name(short_name, index);
        result.push_back(find_result{*p, text, errors, short_name, *p - start + 1});
        if (max_results && loci.size() > max_results) {
          return false;
        }
      }
      for (; q != loci.end() && *q == *p; ++q) {
      }
      p = q;
    }

    return true;
  }

  std::ostream &dump(std::ostream &os) {
    std::string tmp;
    return os;
  }

};

#endif


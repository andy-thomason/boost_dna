
#ifndef DNA_SUFFIX_ARRAY_INCLUDED
#define DNA_SUFFIX_ARRAY_INCLUDED

#include "basic_dna.h"

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
#include <strstream>
#include <algorithm>
#include <numeric>

class dna_suffix_array {
  basic_dna &dna;

  std::vector<uint32_t> index;
  //std::vector<uint8_t> remainder;
  std::vector<uint32_t> addresses;

  template <class _Fn> static void par_for(int limit, _Fn f) {
    //std::vector<std::thread> threads(std::thread::hardware_concurrency());
    std::vector<std::thread> threads(1);

    std::atomic<int> sequence = 0;
    for (auto &t : threads) {
      t = std::thread([&](){
        int seq;
        while ((seq = sequence++) < limit) {
          f(seq);
        }
      });
    }

    for (auto &t : threads) t.join();
  }

  /*template <class... args_t> void log(std::string &str, args_t... args) {
    str.resize(sprintf(nullptr, args...));
    sprintf(str.data(), args...);
    std::cout << str;
  }*/

public:
  dna_suffix_array(basic_dna &dna, int lg_index_size) : dna(dna) {
    index.resize((size_t)1 << lg_index_size);

    const int lg_num_seqs = 4;
    const int num_seqs = 1 << lg_num_seqs;
    size_t counts[num_seqs];
    std::fill(std::begin(counts), std::end(counts), 0);
    std::fill(std::begin(index), std::end(index), 0xffffffff);

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
      //printf("seq=%x num_values=%10d\n", seq, (int)sorter.size());

      size_t addr_idx = seq == 0 ? (size_t)0 : starts[seq-1];
      size_t idx_idx = (size_t)seq << (lg_index_size - lg_num_seqs);
      size_t last_idx_idx = (size_t)(seq+1) << (lg_index_size - lg_num_seqs);
      for (auto p : sorter) {
        while ((p.first >> (64-lg_index_size)) >= idx_idx) {
          //std::cout << "idx " << idx_idx << " " << addr_idx << "\n";
          index[idx_idx++] = (uint32_t)addr_idx;
        }
          //std::cout << "addr_idx " << addr_idx << " " << p.second << "\n";
        addresses[addr_idx++] = (uint32_t)p.second;
      }
      assert(addr_idx == starts[seq]);
      while (idx_idx < last_idx_idx) {
          //std::cout << "+idx " << idx_idx << " " << addr_idx << "\n";
        index[idx_idx++] = (uint32_t)addr_idx;
      }
    });

    assert(std::find(index.begin(), index.end(), 0xffffffff) == index.end());
  }

  void operator=(dna_suffix_array &rhs) = delete;
  dna_suffix_array(basic_dna &rhs) = delete;



  std::ostream &dump(std::ostream &os) {
    std::string tmp;
    return os;
  }

};

#endif


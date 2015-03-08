
#include <cstdint>
#include <fstream>

#include "../include/suffix_array.h"
#include "../include/dna_database.h"
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <iostream>
#include <sstream>



static std::string &to_dna(std::string &str, uint64_t seed, int seed_size=32) {
  str.resize(seed_size);
  for (int i = 0; i != seed_size; ++i) {
    str[i] = "ACGT"[(seed >> 62)&3];
    seed <<= 2;
  }
  return str;
}

#if 0
bool test_ref() {
  static const char rt_fa_test[] =
    ">22 dna:chromosome chromosome:GRCh38:22:1:50818468:1 REF\n"
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAATTCTTGTGTTTATATAA\n"
    "TAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAA\n"
    "GCATTATAAATACAGAGTGCTAAGTTACTTCACTGTGAAATGTAGTCATATAAAGAACAT\n"
    "ANTAATTATACTGGATTATTTTTAAATGGGCTGTCTAACATTATATTAAAAGG\n"
    ">21 dna:chromosome chromosome:GRCh38:22:1:50818468:1 REF\n"
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAATTCTTGTGTTTATATAA\n"
    "GCATTATAAATACAGAGTGCTAAGTTACTTCACTGTGAAATGTAGTCATATAAAGAACAT\n"
    "TAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAA\n"
    "ANTAATTATACTGGATTATTTTTAAATGGGCTGTCTAACATTATATTAAAAGG\n"
  ;

  std::cout << "1" << std::endl;

  dna_database<> dna;
  parser parser;
  parser.add_fasta(&dna, rt_fa_test, rt_fa_test + sizeof(rt_fa_test)-1);

  /*std::ostringstream ss;
  dna.dump(ss);
  std::cout << ss.str();

  if (ss.str() != std::string(rt_fa_test, rt_fa_test + sizeof(rt_fa_test)-1)) {
    throw(std::exception("dump not the same as input\n"));
  }

  std::string buf;
  static const char test_str[] = "TAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAAGCATTATAAATACAGAGTGCTAAGTTACTTCACTGTGAAATGTAGTCATATAAAGAACAT";
  ss.clear();
  for (auto p = dna.begin(); p != dna.end(); ++p) {
    p->for_each_forward_seed(
      [&ss, &buf](uint64_t offset, uint64_t seed) {
        if (offset >= 120 && offset < 200-31) {
          std::string tstr(test_str+offset-120, test_str+offset-120+32);
          std::cout << std::hex << offset << " " << to_dna(buf, seed) << " " << tstr << "\n";
          if (to_dna(buf, seed) != tstr) {
            throw(std::exception("for_each_forward_seed_group not the same as input\n"));
          }
        }
      }, 32
    );
  }*/

  suffix_array<uint32_t, uint32_t> suf(dna, 6);

  std::vector<find_result> result;
  suf.find(dna, result, "TAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAA", 3, 0);
  for (auto &r : result) {
    std::cout << r.locus_ << " " << r.errors << " " << r.dna << "\n";
  }

  suf.find(dna, result, "GAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAA", 3, 0);
  for (auto &r : result) {
    std::cout << r.locus_ << " " << r.errors << " " << r.dna << "\n";
  }

  suf.find(dna, result, "TTT", 1, 0);
  for (auto &r : result) {
    std::cout << r.locus_ << " " << r.errors << " " << r.dna << "\n";
  }

  return true;
}

/*class test_reads {
public:
  test_reads() {
    static const char rt_fq_test[] =
      "@HWI-13CZTJ:55:AXJ7BCZ:1:1101:10685:4594 1:N:0:TTAGGC\n"
      "CGGATGTCAGAGGGGTGCCTTGGGTAACCTCTGGGACTCAGAAGTGAAAGGGGGCTATTCCTAGTTTTATTGCTATAGCCATTATGATTATTAATGATG\n"
      "+\n"
      "CCCFFFFFHHHHHJJCFHJJJIJJDGIJJJJJIJJIJJJJJIJJCHHIJJJIJHFFFEEEEDDECDEDDDDDDDEDEEDDDDDDEDDEEEDEEEEDDDD\n"
      "@HWI-13CZTJ:55:AXJ7BCZ:1:1101:10685:4594 2:N:0:TTAGGC\n"
      "ACTGAAGCTCGTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTCTAATAGCTATCCTCTTCAACAATATACTCTC\n"
      "+\n"
      "@CCFFFFFHHHHHJJJJJJJJJJJJJIJIJJJJ?HIIJJIJJJJJJJJJJGHIIJJIJHHHHHFFFFFFEEEEEEDDDEDDDEDDDDEDDDDDDDEEEDDD\n"
    ;

    reads<def_read_traits> test(rt_fq_test, rt_fq_test + sizeof(rt_fq_test)-1);

    std::ostringstream ss;
    test.dump(ss);
    std::cout << ss.str();

    if (ss.str() != std::string(rt_fq_test, rt_fq_test + sizeof(rt_fq_test)-1)) {
      throw(std::exception("dump not the same as input\n"));
    }

    std::string buf;
    static const char *test_str[] = {
      "CGGATGTCAGAGGGGTGCCTTGGGTAACCTCTGGGACTCAGAAGTGAAAGGGGGCTATTCCTAGTTTTATTGCTATAGCCATTATGATTATTAATGATG",
      "AAAAAAGCTCGTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTCTAATAGCTATCCTCTTCAACAATATACTCTC",
    };

    for (int seed_size = 64; seed_size <= 32; seed_size <<= 1) {
      for (auto item : test) {
        item.for_each_forward_seed(
          [&ss, &buf, seed_size](uint64_t addr, uint64_t seed) {
            uint64_t idx = addr >= 99 ? addr - 99 : addr;
            const char *p = test_str[addr >= 99];
            std::string tstr(p + idx, p + idx + seed_size);
            std::cout << std::hex << addr << " " << to_dna(buf, seed, seed_size) << " " << tstr << "\n";
            if (to_dna(buf, seed, seed_size) != tstr) {
              throw(std::exception("for_each_forward_seed_group not the same as input\n"));
            }
          }, seed_size
        );
      }
    }
    for (auto item : test) {
      item.for_each_forward_seed(
        [](uint64_t addr, uint64_t seed) {
          volatile uint64_t a = addr, s = seed;
        }, 32
      );
    }
  }
};*/


bool test_file() {
  using namespace boost::interprocess;
  file_mapping fa_file("C:/projects/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", read_only);

  //file_mapping fa_file("C:/projects/test/chr21_p1.fa", read_only);
  mapped_region region(fa_file, read_only, 0, 0);

  char *begin = (char*)region.get_address();
  char *end = begin + region.get_size();
  dna_database<> dna;
  parser p;
  //p.add_fasta(&dna, begin, end);

  if (0) {
    std::ofstream wfile("dna.dat", std::ios_base::binary);
    dna.write(wfile);
  }

  {
    std::ifstream rfile("dna.dat", std::ios_base::binary);
    dna.read(rfile);
  }

  if (0) {
    suffix_array<uint32_t, uint32_t> suf;
    if (0) {
      suf = suffix_array<uint32_t, uint32_t>(dna, 24);
      std::ofstream wfile("suf.dat", std::ios_base::binary);
      suf.write(wfile);
      return true;
    }
    if (0) {
      std::ifstream rfile("suf.dat", std::ios_base::binary);
      suf.read(rfile);
    }

    std::vector<find_result> result;
    suf.find(dna, result, "GACATCATCCTGTACGCGTC", 1, 0);
    for (auto &r : result) {
      std::cout << r.chromosome << " " << r.offset << " " << r.errors << " " << r.dna << "\n";
    }
    return true;
  }
}

#endif

void make_from_file(const char *input, const char *dna_file, const char *suffix_file) {
  using namespace boost::interprocess;
  file_mapping fa_file(input, read_only);
  mapped_region region(fa_file, read_only, 0, 0);

  char *begin = (char*)region.get_address();
  char *end = begin + region.get_size();
  parser p;
  dna_database_def dna;
  p.add_fasta(&dna, begin, end);

  dna.write(std::ofstream(dna_file, std::ios_base::binary));

  suffix_array_def suf(dna, 24);
  suf.write(std::ofstream(suffix_file, std::ios_base::binary));
}

void test_one(const char *dna_file, const char *suffix_file, const char *sequence, int max_distance) {
  dna_database_def dna;
  dna.read(std::ifstream(dna_file, std::ios_base::binary));

  std::vector<find_result> results;
  long long t1 = __rdtsc();
  //dna.find(results, sequence, max_distance, 0);
  long long t2 = __rdtsc();
  for (auto &r : results) {
    std::cout << r << "\n";
  }

  boost::interprocess::file_mapping mapping(suffix_file, boost::interprocess::read_only);
  boost::interprocess::mapped_region region(mapping, boost::interprocess::read_only, 0, 0);
  const char *begin = (const char*)region.get_address();
  const char *end = begin + region.get_size();

  suffix_array_mapped suf;
  suf.map(begin, end);

  if (0) {
    std::cout << "verifying\n";
    bool res = suf.verify(dna);
    std::cout << "verified" << res << "\n";
  }

  std::cout << t2-t1 << "\n";

  long long t3 = __rdtsc();
  suf.find(dna, results, sequence, max_distance, 0);
  long long t4 = __rdtsc();
  for (auto &r : results) {
    std::cout << r << "\n";
  }

  std::cout << t4-t3 << "\n";
}

class allocator_base {
public:
};

#define noexcept
template <class T> class allocator {
public:
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T value_type;
  template <class U> struct rebind { typedef allocator<U> other; };

  allocator(allocator_base *base = nullptr) noexcept : base(base) {}
  allocator(const allocator&x) noexcept { base = x.base; }
  template <class U> allocator(const allocator<U> &x) noexcept { base = x.base; }
  ~allocator() {}
  pointer address(reference x) const noexcept { return &x; }
  const_pointer address(const_reference x) const noexcept { return &x; }
  pointer allocate(size_type count, allocator<void>::const_pointer hint = 0) { return (pointer)malloc(count * sizeof(T)); }
  void deallocate(pointer p, size_type n) noexcept { free((void*)p); }
  size_type max_size() const noexcept { return 0xffffffff; }
  template<class U, class... Args> void construct(U* p, Args&&... args) { new ((void*)p)U(args...); }
  template <class U> void destroy(U* p) { p->~U(); }

  allocator_base *base;
};
#undef noexcept

template <> struct allocator<void> {
  typedef void* pointer;
  typedef const void* const_pointer;
};

void test_vector() {
  //std::vector<uint64_t> allocated(256);
  //allocator<int>{}
  //std::vector<int, allocator<int>> mv();



}

template <class _FwdIter, class _OutTy> void make_suffix_array(_FwdIter begin, _FwdIter end, _OutTy dest) {
  typedef decltype(*begin) value_type;
  typedef decltype(*dest) index_type;
  std::iota(dest, dest + ((end - begin) + 1), 0);
}

void test_suffix_array() {
  std::string test = "abracadabra";
  std::vector<std::uint8_t> sa(test.size() + 1);
  make_suffix_array(test.begin(), test.end(), sa.begin());
}


int main() {

  try {
    test_suffix_array();
    //make_from_file("C:/projects/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", "GRCh38.dna", "GRCh38.suf");
    //make_from_file("C:/projects/test/chr21_p1.fa", "chr21.dna", "chr21.suf");

    //test_one("chr21.dna", "chr21.suf", "CCAATGCCTAGGGAGATTTCTAGGTCCTCTGTTCCTTGCTGACCTCCAAT", 3);
    //test_one("GRCh38.dna", "GRCh38.suf", "CCAATGCCTAGGGAGATTTCTAGGTCCTCTGTTCCTTGCTGACCTCCAAT", 3);
    //test_file();
    //test_ref();
    std::cout << "Passed:\n";
  } catch(std::exception e) {
    std::cout << "Failed: " << e.what() << "\n";
    return 1;
  }
}

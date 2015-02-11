
#include <cstdint>

struct hex {
  char str[17];

  template <class _Ty> hex(_Ty n) {
    uint64_t value = (uint64_t)n << (64 - sizeof(n)*8);
    for (int i = 0; i != sizeof(n)*2; ++i) {
      str[i] = "0123456789abcdef"[value >> 60];
      value <<= 4;
    }
    str[sizeof(n)*2] = 0;
  }

  operator const char *() { return str; }
};

#include "../include/suffix_array.h"
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

  dna_database dna;
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

  suffix_array suf(dna, 6);

  std::vector<suffix_array::find_result> result;
  suf.find(result, "TAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAA", 3, 0);
  for (auto &r : result) {
    std::cout << r.locus_ << " " << r.errors << " " << r.dna << "\n";
  }

  suf.find(result, "GAAGATGTCCTATAATTTCTGTTTGGAATATAAAATCAGCAACTAATATGTATTTTCAAA", 3, 0);
  for (auto &r : result) {
    std::cout << r.locus_ << " " << r.errors << " " << r.dna << "\n";
  }

  suf.find(result, "TTT", 1, 0);
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
  //file_mapping fa_file("C:/projects/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", read_only);
  file_mapping fa_file("C:/projects/test/chr21_p1.fa", read_only);

  mapped_region region(fa_file, read_only, 0, 0);
  char *begin = (char*)region.get_address();
  char *end = begin + region.get_size();
  dna_database dna;
  parser p;
  p.add_fasta(&dna, begin, end);
  return true;
}

int main() {
  try {
    //test_file();
    test_ref();
    std::cout << "Passed:\n";
  } catch(std::exception e) {
    std::cout << "Failed: " << e.what() << "\n";
    return 1;
  }
}

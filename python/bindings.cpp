
#include <cstdint>

#include <boost/python.hpp> //relative to ENV system import path
#include <boost/python/list.hpp> // Header file for PyList
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include "../include/suffix_array.h" //means relative path import

bool build(boost::python::list input_files, char* output_file)
{
  using namespace boost::interprocess;
  //file_mapping fa_file("C:/projects/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", read_only);

  printf("%d", boost::python::len(input_files));

  /*file_mapping fa_file("C:/projects/test/chr21_p1.fa", read_only);

  mapped_region region(fa_file, read_only, 0, 0);
  char *begin = (char*)region.get_address();
  char *end = begin + region.get_size();
  dna_database dna;
  parser p;
  p.add_fasta(&dna, begin, end);
  */
  return true;
}

void find_inexact()
{
    BOOST_PYTHON_MODULE(suffix_array)
    using namespace boost::python;
    def("build", build);
    def("find_inexact", find_inexact);
}



#include <cstdint>
#include <fstream>
#include <future>

#define HAVE_ROUND
#include <boost/python.hpp> //relative to ENV system import path
#include <boost/python/list.hpp> // Header file for PyList
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include "../include/dna_database.h" //means relative path import
#include "../include/suffix_array.h" //means relative path import

int INDEX_POWER = 24;

bool
build(boost::python::list input_files, char* output_file)
{
  using namespace boost::interprocess;
  using namespace boost::python;

  //file_mapping fa_file("C:/projects/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", read_only);
  size_t num_files = len(input_files);

  dna_database_def dna;

  for (size_t i=0; i < num_files; ++i ) {
    char* path_to_file = extract<char*>(input_files[i]);
    std::cerr << path_to_file << "\n";
    //add to the db
    file_mapping fa_file(path_to_file, read_only); //open file
    mapped_region region(fa_file, read_only, 0, 0); //const map window
    char *begin = (char*)region.get_address(); //set start of window
    char *end = begin + region.get_size(); //set end of window 0 = size of file
    //init DNA DB

    dna.begin_section();
    std::string chr_name = path_to_file;
    // get a name for the molecule
    dna.add_name(chr_name.data(), chr_name.data()+chr_name.size()); //don't try this at home
    dna.add_dna(begin, end); //add to packed db all bp to ints
    dna.end_section();  //reads entire file
  }

  // Now need to construct the SA (and do the big sort)
  suffix_array_def sa(dna, INDEX_POWER); //pointer to the db and resolution of index

  // dump it to a file
  std::ofstream os(output_file, std::ios::binary);
  sa.write(os);
  dna.write(os);
  os.close();

  return output_file;
}

boost::python::list
find_inexact(char* cpp_sa, char* pattern, int maxmismatch, int limit)
{
  suffix_array_def sa;
  dna_database_def dna;
  std::ifstream is(cpp_sa, std::ios::binary);
  sa.read(is);
  dna.read(is);
  is.close();

  using namespace boost::python;
  list results;
  // find_result is the element of the result vector
  //  it has our data in it
  std::vector<find_result> result;

  sa.find(dna, result, pattern,  maxmismatch, limit);
  // copy
  for(auto r : result) //iterate over C++ vector and make a pyobject
  {
    // need to make a pyobject to send back to Python
    results.append(object(r.chromosome()));
  }
  return results;
}

BOOST_PYTHON_MODULE(suffix_array)
{
  using namespace boost::python;
  def("build", build);
  def("find_inexact", find_inexact);
}


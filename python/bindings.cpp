
#include <cstdint>

#include <boost/python.hpp> //relative to ENV system import path
#include <boost/python/list.hpp> // Header file for PyList
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include "../include/suffix_array.h" //means relative path import

int INDEX_POWER = 24;

dna_database dna; //naughty so fix later
suffix_array* sa; // our suffix array

bool build(boost::python::list input_files, char* output_file)
{
  using namespace boost::interprocess;
  using namespace boost::python;

  uint8_t translate[256];
  std::fill(std::begin(translate), std::end(translate), 4);
  translate['\r'] = 8;
  translate['\n'] = 8;
  translate['A'] = 0;
  translate['C'] = 1;
  translate['G'] = 2;
  translate['T'] = 3;

  //file_mapping fa_file("C:/projects/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", read_only);
  int num_files = len(input_files);

  for(int i=0;i<num_files;++i)
  {
    char* path_to_file = extract<char*>(input_files[i]);
    std::cerr<<path_to_file<<"\n";
    //add to the db
    file_mapping fa_file(path_to_file, read_only); //open file
    mapped_region region(fa_file, read_only, 0, 0); //const map window
    char *begin = (char*)region.get_address(); //set start of window
    char *end = begin + region.get_size(); //set end of window 0 = size of file
    //init DNA DB
    // parser p; not required if the files are flat
    //p.add_fasta(&dna, begin, end);
    dna.begin_section();
    std::string chr_name = path_to_file;
    // get a name for the molecule
    dna.add_name(chr_name.data(), chr_name.data()+chr_name.size()); //don't try this at home
    dna.add_dna(begin, end, translate); //add to packed db all bp to ints
    dna.end_section();  //reads entire file
  }
  // Now need to construct the SA (and do the big sort)
  sa = new suffix_array(dna, INDEX_POWER); //pointer to the db and resolution of index


  //dump it to a file
  // TODO add a load from disk and write to disk

  return output_file;
}

boost::python::list
find_inexact(char* cpp_sa, char* pattern, int maxmismatch, int limit)
{

  using namespace boost::python;
  list results;
  //find_result is the element of the result vector
  // it has our data in it
  std::vector<suffix_array::find_result> result;

  sa->find(result, pattern,  maxmismatch, limit);
  //copy
  for(auto r:result) //iterate over C++ vector and make a pyobject
  {
    //need to make a pyobject to send back to Python
    results.append(object(r.chromosome));
  }
  return results;
}

BOOST_PYTHON_MODULE(suffix_array)
{
    using namespace boost::python;
    def("build", build);
    def("find_inexact", find_inexact);
}


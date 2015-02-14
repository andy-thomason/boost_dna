"""
This module tests the inexact match search for the boost dna module
It will assume that a genome is present and installed.

SCOPE: functional
"""

import pytest
import os

import suffix_array  # This module should be generated by Boost


PATH_TO_GENOME = '/home/rileyd/biodata/genomes/Human/GRCh38'



#@pytest.fixture
#def cpp_sa(request):
def cpp_sa():
    """
    Construct a new suffix array and let it hang around for the test session
    :param scope:
    :return:
    """
    # The files of our genome
    sequence_file_paths = []

    # the path to our index
    index_path = os.path.join(PATH_TO_GENOME, 'index/suffix_array.sa')
    for filename in os.listdir(PATH_TO_GENOME):
        if filename.endswith('.flat'):  # edit to change to flat files
            sequence_uri = os.path.join(PATH_TO_GENOME, filename)
            sequence_file_paths.append(sequence_uri)

    # invoke the build code
    error = suffix_array.build(sequence_file_paths, index_path)

    def clean_up():  # REMOVE THE INDEX EACH RUN
        os.remove(index_path)

    #request.addfinalizer(clean_up)

    return index_path


# put some tests here
def test_simple_case(cpp_sa):
    """
    Test a relatively unique pattern for a simple basic case
    :return:
    """
    genome_version = "GRCh38.p2"
    pattern = "GACATCATCCTGTACGCGTC"
    chromosome_name = "chr5"
    offset_in_chromosome = 177093184 - 17

    offtarget_one = ("GGCATCACCCTGTACCCGTC", "chr5", 139828809 - 17)
    offtarget_two = ("GACATCATCCTGGAAGGGTC", "chr1", 29314933 - 17)
    maxmismatch = 3
    limit = 0

    results = suffix_array.find_inexact(cpp_sa, pattern, maxmismatch, limit)
    # Does not filter out "self" hit
    assert len(results) == 3

if __name__ == '__main__':
    print "Starting test"
    cpp_sa()





#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Boost DNA Python Wrapper. Used for searching the genome using a suffix array.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import functools
import logging
import itertools
import operator

import suffix_array


def build(sequence_file_paths, index_path="genome.sa", ):
    """
    Driver to build a Suffix Array index over a collection of DNA sequence files
    :return:
    """
    error = suffix_array.build(sequence_file_paths, index_path)
    if error:
        raise Exception()
    return False


def find_inexact(sa_path, pattern, maxmismatch=4, limit=0):
    """Find up to limit occurrences of a pattern with up to maxmismatch
    mismatches
    """
    results = suffix_array.find_inexact(sa_path, pattern, maxmismatch, limit)



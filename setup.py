from setuptools import setup, Extension

import os

os.environ["CXX"] = "clang++"

setup(
    name='suffix_array',
    version='0.0.1',
    install_requires=[
        'pytest',
    ],
    ext_modules=[
        Extension('suffix_array',
                  ['python/bindings.cpp'],
                  include_dirs=['include', '/usr/include/python2.7'],
                  library_dirs=['/usr/lib/x86_64-linux-gnu/'],
                  libraries=['boost_python'],
                  extra_compile_args=['-std=c++11']
        )
    ],
)

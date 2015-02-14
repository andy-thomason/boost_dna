cmake_minimum_required(VERSION 2.8.4)
project(boost_dna)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "DEBUG")
  #SET(CMAKE_BUILD_TYPE "RELEASE")
  #SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
  #SET(CMAKE_BUILD_TYPE "MINSIZEREL")
ENDIF()

FIND_PACKAGE(Boost 1.45.0)
FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(Boost COMPONENTS python)

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}"
    "/usr/include/python2.7"
    "include"
  )
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.45.0 COMPONENTS python)

  ADD_LIBRARY(suffix_array SHARED python/bindings.cpp)
  TARGET_LINK_LIBRARIES(suffix_array ${Boost_LIBRARIES})
  set_target_properties(suffix_array PROPERTIES PREFIX "" )
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
  ADD_DEFINITIONS("-Wall")
ELSE()
  MESSAGE(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
ENDIF()





#clang++ -O2 -std=c++11 -fPIC -shared -L/usr/lib/x86_64-linux-gnu/ -lboost_python-py27 -I /usr/include/python2.7 python/bindings.cpp -o suffix_array.so

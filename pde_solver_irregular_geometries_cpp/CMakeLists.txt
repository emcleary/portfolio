cmake_minimum_required(VERSION 3.20)
project(Solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -fopenmp -Wno-array-bounds")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "-g -O0 -Wall -Wno-array-bounds")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(src)

find_package(gsl)
find_package(cblas)

add_subdirectory(src)
add_subdirectory(examples)

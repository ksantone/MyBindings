cmake_minimum_required(VERSION 3.4)
project(denovo_sequencing LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 11)

add_subdirectory(extern/pybind11)

pybind11_add_module(denovo_sequencing DeNovoSequencingAlgorithm.cpp DeNovoSequencingAlgorithm.h)
target_link_libraries(denovo_sequencing PUBLIC pybind11::embed)
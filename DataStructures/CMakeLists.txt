project(data_structures_project)

add_library(data_structures AminoAcids.cpp AminoAcids.h FragmentWeightMatrix.cpp Peak.cpp SpectrumGraph.h SpectrumGraph.cpp bindings/)

include_directories(../extern/pybind11)
target_link_libraries(data_structures PUBLIC pybind11::embed)
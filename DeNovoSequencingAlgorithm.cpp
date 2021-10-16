#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/eval.h>
#include <pybind11/iostream.h>
#include "DataStructures/AminoAcids.h"
#include "DataStructures/FragmentWeightMatrix.h"
#include "DataStructures/SpectrumGraph.h"
#include "DeNovoSequencingAlgorithm.h"

namespace py = pybind11;

std::mutex my_mutex;

void print_func(int id) {
    my_mutex.lock();
    std::cout << "Printing from thread: " << std::to_string(id) << "\n";
    my_mutex.unlock();
}

void DeNovoSequencingAlgorithm::run_denovo() {
    std::string spectraFileName = "/Users/kassimsantone/Desktop/spectrum_list.txt";
    std::ifstream spectraFile(spectraFileName);
    int first = 0;
    std::list<SpectrumGraph> spectrumGraphs, filteredSpectrumGraphs;
    std::string spectralData;
    while (getline(spectraFile, spectralData) && first <= 1) {
        if (first == 1) {
            spectrumGraphs.push_back(SpectrumGraph(spectralData));
        }
        first += 1;
    }
    std::cout << "Number of spectra is: " << first << std::endl;
    std::list<SpectrumGraph>::iterator it;
    std::list<FragmentWeightMatrix> fragmentWeightMatrices;
    for (it = spectrumGraphs.begin(); it != spectrumGraphs.end(); it++) {
        fragmentWeightMatrices.push_back(FragmentWeightMatrix(*it, it->spectralPeaks.size()));
    }
}

int main() {
    DeNovoSequencingAlgorithm deNovo;
    deNovo.DeNovoSequencingAlgorithm::run_denovo();
    return 0;
}
/*    std::thread t0(print_func, 0);
    std::thread t1(print_func, 1);
    std::thread t2(print_func, 2);
    std::thread t3(print_func, 3);

    t0.join();
    t1.join();
    t2.join();
    t3.join();
    return 0;
    std::vector<std::thread> threads;
    int sum = 0;
    for (int j = 0; j < 48; j++) {
        threads.push_back(std::thread([&]() { sum += add(2, 4); }));
    }
    for (int j = 0; j < 48; j++) {
        threads.at(j).join();
    }
    return sum;
}*/

PYBIND11_MODULE(denovo_sequencing, handle) {
	handle.doc() = "This is the DeNovo sequencing algorithm.";
	handle.def("main", 
    []() {
        pybind11::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
        main();
    });

	/*py::class_<FragmentWeightMatrix>(
		handle, "FragmentWeightMatrix"
		);*/
}
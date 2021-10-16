#include "FragmentWeightMatrix.h"
#include <thread>
#include <condition_variable>
#include <deque>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <mutex>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/eval.h>
#include <pybind11/iostream.h>

std::mutex mtx;

std::condition_variable current_cv;
std::condition_variable threads_completing_cv;
std::vector<bool> triangle_complete;
std::mutex m1, m2, m3;
int threads_executing = 0;
int rows_completed = 0;

void FragmentWeightMatrix::traverse_fragment_matrix(int id) {
	std::cout << "Traversing fragment matrix: " << id << std::endl;
	threads_executing -= 1;
	std::cout << "Threads executing after traversal: " << threads_executing << std::endl;
}

FragmentWeightMatrix::FragmentWeightMatrix(const FragmentWeightMatrix &fragmentWeightMatrix) {
	std::cout << "Copied!" << std::endl;
	FragmentWeightMatrix::spectrumPeaks = fragmentWeightMatrix.spectrumPeaks;
	FragmentWeightMatrix::edges = fragmentWeightMatrix.edges;
	std::cout << "Done copy..." << std::endl;
}

FragmentWeightMatrix::FragmentWeightMatrix(SpectrumGraph spectrumGraph, int spectrumPeaks) {
	FragmentWeightMatrix::spectrumPeaks = spectrumPeaks;
	edges = spectrumGraph.spectralEdges;
	std::cout << "Number of peaks is: " << spectrumPeaks << std::endl;
	for (int i = 0; i < spectrumPeaks * spectrumPeaks; i++) {
		fragmentWeightMatrix.push_back(0);
	}
	for (int i = 1; i < spectrumPeaks; i++) {
		fragmentWeightMatrix.at(i + i * spectrumPeaks) = -std::numeric_limits<double>::infinity();
	}
	std::cout << "Just before thread execution!" << std::endl;
	std::thread traverse_thread;
	for (int i = 0; i <= spectrumPeaks; i++) {
		std::cout << "On thread number: " << i+1 << "." << std::endl;
		traverse_thread = std::thread([&]() {this->FragmentWeightMatrix::traverse_fragment_matrix(i);});
		traverse_thread.join();
		std::cout << "Thread number " << i+1 << " completed." << std::endl;
		traverse_thread.detach();
	}
}
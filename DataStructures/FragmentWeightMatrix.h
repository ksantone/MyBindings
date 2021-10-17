#include "SpectrumGraph.h"
#include <array>
#include <condition_variable>
#include <mutex>

#ifndef FragmentWeightMatrix_HEADER
#define FragmentWeightMatrix_HEADER

class FragmentWeightMatrix {
public:
	std::condition_variable current_cv;
	std::condition_variable threads_completing_cv;
	std::vector<bool> triangle_complete;
	std::mutex m1, m2, m3;
	int threads_executing = 0;
	int rows_completed = 0;
	int spectrumPeaks;
	int num_of_threads = 4;
	int smallest_available_thread = 8;
	std::vector<bool> row_complete;
	std::vector<int> available_indices;
	std::vector<double> fragmentWeightMatrix;
	void construct_triangle(int);
	void compute_row(int);
	void traverse_fragment_matrix(int);
	double findMaximumSubfragment(std::vector<double>, int, int, std::map<std::pair<int, int>, std::string>, int);
	std::vector<std::string> glue_edge_to_all_sequences(std::string, std::vector<std::string>, bool);
	std::vector<std::string> find_sequence_from_fragment_position(std::vector<double>, std::map<std::pair<int, int>, std::string>, std::vector<std::pair<int, int> >, double, std::pair<int, int>, int, std::string);
	
	std::map<std::pair<int, int>, std::string> edges;
	FragmentWeightMatrix(SpectrumGraph, int);
	FragmentWeightMatrix(const FragmentWeightMatrix &FragmentWeightMatrix);
	
};

#endif
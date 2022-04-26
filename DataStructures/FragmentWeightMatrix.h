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
	std::vector<std::string> glue_edge_to_all_sequences(std::string, std::pair<int, int>, std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > >, bool);
	std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > find_b_sequence_from_fragment_position(std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > >, std::vector<double>, std::map<std::pair<int, int>, std::string>, std::vector<std::pair<int, int> >, int, std::vector<std::pair<int, int> >, std::string, std::vector<std::vector<std::pair<std::string, int> > >, std::vector<Peak>);
	std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > find_y_sequence_from_fragment_position(std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > >, std::vector<double>, std::map<std::pair<int, int>, std::string>, std::vector<std::pair<int, int> >, int, std::vector<std::pair<int, int> >, std::string, std::vector<std::vector<std::pair<std::string, int> > >, std::vector<Peak>);
	std::vector<std::vector<std::pair<int, int> > > glue_edge_to_position_vector(std::pair<int, int>, std::vector<std::vector<std::pair<int, int> > >, bool);
	std::map<std::pair<int, int>, std::string> edges;
	FragmentWeightMatrix(SpectrumGraph, int, std::map<std::pair<int, int>, std::string>);
	FragmentWeightMatrix(const FragmentWeightMatrix &FragmentWeightMatrix);
	
};

#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/eval.h>
#include <pybind11/iostream.h>
#include "DataStructures/FragmentWeightMatrix.h"
#include "DataStructures/SpectrumGraph.h"
#include "DeNovoSequencingAlgorithm.h"

namespace py = pybind11;

// Global variables and classes

std::mutex my_mutex;
std::map <std::string, float> amino_acids = { {"A", 71.03711}, {"R", 156.10111}, {"N", 114.04293}, {"D", 115.02694}, {"C", 103.00919}, {"E", 129.04259}, {"Q", 128.05858 },
        {"G", 57.02146}, {"H", 137.05891}, {"I/L", 113.08406}, {"K", 128.09496}, {"M", 131.04049}, {"F", 147.06841}, {"P", 97.05276},
        {"S", 87.03203}, {"T", 101.04768}, {"W", 186.07931}, {"Y", 163.06333}, {"V", 99.06841} };
std::mutex mtx;
std::condition_variable current_cv;
std::condition_variable threads_completing_cv;
std::vector<bool> triangle_complete;
std::mutex m1, m2, m3;
int threads_executing = 0;
int rows_completed = 0;

class AminoAcids {
public:
    static std::map <std::string, float> amino_acids;
    std::map<std::string, float> getAminoAcids();
    float getAminoAcidMass(std::string amino_acid);
    AminoAcids();
};

AminoAcids::AminoAcids() {

}

float round_var(float var)
{
	float value = (int)(var * 10 + .5);
	return (float)value / 10;
}

// SpectrumGraph-related methods

std::map<std::pair<int, int>, std::string> SpectrumGraph::computeSpectralEdges(std::vector<Peak> spectralPeaks) {
	std::map<std::pair<int, int>, std::string> output_edges;
	for (int i = 0; i < spectralPeaks.size(); i++) {
		// Check if spectralPeaks differences are in dictionary (maybe make dictionary round down to near 1/10^n then check if either equal
		// or less than actual value by 1/10^n
		for (int j = i + 1; j < spectralPeaks.size(); j++) {
			float peak_diff = spectralPeaks.at(j).mz - spectralPeaks.at(i).mz;
			if (peak_diff >= 187) {
				break;
			}
			else if (peak_diff <= 56) {
				continue;
			}
			else {
				std::map<std::string, float>::iterator it;
				for (it = amino_acids.begin(); it != amino_acids.end(); it++) {
					if (std::abs((it->second)-peak_diff)<0.1) {
						//output_edges.emplace_back(std::pair<Peak, Peak>(spectralPeaks.at(i), spectralPeaks.at(j)));
						output_edges[std::pair <int, int>(i, j)] = it->first;
						// Sometimes there are multiple amino acids that match (unless we reduce to say 0.4, but for now let's just
						// take the first one)
						break;
					}
				}
			}
		}
	}
	return output_edges;
}

SpectrumGraph::SpectrumGraph(std::string spectralData) {
	std::istringstream iss(spectralData);
	std::string del = ",";
	std::string peak_info;
	std::stringstream ss;
	float mz, rt, intensity;
	iss >> rt;
	while (iss >> peak_info) {
		int start = 0;
		int end = peak_info.find(del);
		ss << peak_info.substr(start, end - start) << std::endl;
		mz = std::stof(ss.str());
		start = end + del.size();
		end = peak_info.size();
		ss.str("");
		ss.clear();
		ss << peak_info.substr(start, end - start) << std::endl;
		intensity = std::stof(ss.str());
		spectralPeaks.push_back(Peak(mz, intensity, rt));
		ss.str("");
		ss.clear();
	}
	std::vector<Peak>::const_iterator firstPeak = spectralPeaks.begin();
	std::vector<Peak>::const_iterator middlePeak = spectralPeaks.begin() + (int)(spectralPeaks.size() / 2);
	std::vector<Peak> newSpectralPeaks = std::vector<Peak>(firstPeak, middlePeak);
	spectralEdges = computeSpectralEdges(spectralPeaks);
	spectralPeaks = newSpectralPeaks;
}

Peak::Peak(float mz, float intensity, float rt)
{
	this->mz = mz;
	this->intensity = intensity;
	this->rt = rt;
}

// FragmentWeightMatrix-related methods

double FragmentWeightMatrix::findMaximumSubfragment(std::vector<double> fragmentWeightMatrix, int i, int j, std::map<std::pair<int, int>, std::string> edges, int spectrumPeaks) {
	double currentMaxfragment = -1;
	std::vector<std::pair<int, int> > edges_without_amino_acids;
	transform(edges.begin(), edges.end(), std::back_inserter(edges_without_amino_acids), [&](std::pair<std::pair<int, int>, std::string> edge) { return edge.first; });
	std::pair<int, int> final_pair;
	std::string amino_acid;
	for (int k = j - 1; k >= 0; k--) {
		double current_entry = fragmentWeightMatrix.at(i * spectrumPeaks + k);
		if (current_entry > currentMaxfragment && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, j))) {
			currentMaxfragment = current_entry + 1;
			final_pair = std::pair<int, int>(i, k);
			amino_acid = edges[std::pair<int, int>(k, j)];
		}
	}
	for (int k = i - 1; k >= 0; k--) {
		double current_entry = fragmentWeightMatrix.at(k * spectrumPeaks + j);
		if (current_entry > currentMaxfragment && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, i))) {
			currentMaxfragment = current_entry + 1;
			final_pair = std::pair<int, int>(k, j);
			amino_acid = edges[std::pair<int, int>(k, i)];
		}
	}
	if (currentMaxfragment > -1) {
		mtx.lock();
		std::cout << "Computed fragment at position: (" << i << ", " << j << ") to have " << currentMaxfragment << " peaks." << std::endl;
		mtx.unlock();
	}
	return currentMaxfragment;
}

std::vector<std::string> FragmentWeightMatrix::glue_edge_to_all_sequences(std::string amino_acid, std::vector<std::string> subpeptides, bool left) {
	std::vector<std::string> new_subpeptides;
	for (int i = 0; i < subpeptides.size(); i++) {
		if (left) {
			new_subpeptides.push_back(amino_acid + subpeptides.at(i));
		}
		else {
			new_subpeptides.push_back(subpeptides.at(i) + amino_acid);
		}
	}
	return new_subpeptides;
}

std::vector<std::string> FragmentWeightMatrix::find_sequence_from_fragment_position(std::vector<double> fragmentWeightMatrix, std::map<std::pair<int, int>, std::string> edges, std::vector<std::pair<int, int> > edges_without_amino_acids, double max_num_of_peaks, std::pair<int, int> optimal_fragment_position, int spectrumPeaks, std::string prefix_suffix_connector) {
	std::vector<std::string> empty_string_vector;
	if (max_num_of_peaks == 0) {
		empty_string_vector.push_back(prefix_suffix_connector);
		return empty_string_vector;
	}
	int i = optimal_fragment_position.first, j = optimal_fragment_position.second;
	double max_prev_frag = -std::numeric_limits<double>::infinity();
	std::pair<int, int> prev_frag_pair;
	bool found = false;
	std::vector<std::string> vector_1, vector_2;
	for (int k = i - 1; k >= 0; k--) {
		if (fragmentWeightMatrix.at(k * spectrumPeaks + j) > max_prev_frag && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, i))) {
			max_prev_frag = fragmentWeightMatrix.at(k * spectrumPeaks + j);
			prev_frag_pair = std::pair<int, int>(k, i);
			found = true;
			vector_1 = FragmentWeightMatrix::glue_edge_to_all_sequences(edges[prev_frag_pair], find_sequence_from_fragment_position(fragmentWeightMatrix, edges, edges_without_amino_acids, max_num_of_peaks - 1, std::pair<int, int>(k, j), spectrumPeaks, prefix_suffix_connector), true);
		}
	}
	for (int k = j - 1; k >= 0; k--) {
		if (fragmentWeightMatrix.at(i * spectrumPeaks + k) > max_prev_frag && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, j))) {
			max_prev_frag = fragmentWeightMatrix.at(i * spectrumPeaks + k);
			prev_frag_pair = std::pair<int, int>(k, j);
			found = true;
			vector_2 = FragmentWeightMatrix::glue_edge_to_all_sequences(edges[prev_frag_pair], find_sequence_from_fragment_position(fragmentWeightMatrix, edges, edges_without_amino_acids, max_num_of_peaks - 1, std::pair<int, int>(i, k), spectrumPeaks, prefix_suffix_connector), false);
		}
	}
	if (!found) {
		return empty_string_vector;
	}
	vector_1.insert(vector_1.end(), vector_2.begin(), vector_2.end());
	return vector_1;
}

void FragmentWeightMatrix::compute_row(int row) {
	for (int j = row + 1; j < spectrumPeaks; j++) {
		// Compute on row to right of diagonal
		double temp = FragmentWeightMatrix::findMaximumSubfragment(fragmentWeightMatrix, row, j, edges, spectrumPeaks);
		fragmentWeightMatrix.at(row * spectrumPeaks + j) = temp;
		// Follows that column below diagonal is the same since matrix is symmetric
		fragmentWeightMatrix.at(j * spectrumPeaks + row) = fragmentWeightMatrix.at(row * spectrumPeaks + j);
	}
	threads_executing -= 1;
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
	pybind11::gil_scoped_release release;
	std::thread traverse_thread;
	for (int i = 0; i <= spectrumPeaks; i++) {
		std::cout << "On thread number: " << i+1 << "." << std::endl;
		traverse_thread = std::thread([&]() {this->FragmentWeightMatrix::compute_row(i);});
		traverse_thread.join();
		std::cout << "Thread number " << i+1 << " completed." << std::endl;
	}
	pybind11::gil_scoped_acquire acquire;
	double max_num_of_peaks = 0;
	std::pair<int, int> optimal_fragment_position;
	std::vector<std::pair<int, int> > edges_without_amino_acids;
	transform(edges.begin(), edges.end(), std::back_inserter(edges_without_amino_acids), [&](std::pair<std::pair<int, int>, std::string> edge) { return edge.first; });
	std::vector<std::string> prefix_suffix_connectors;
	std::vector<std::pair<int, int> > optimal_fragment_positions;
	std::map<std::pair<int, int>, double> fragment_to_num_of_peaks;
	for (int i = 0; i < spectrumPeaks; i++) {
		for (int j = 0; j < spectrumPeaks; j++) {
			if (fragmentWeightMatrix.at(i * spectrumPeaks + j) == -std::numeric_limits<double>::infinity()) {
				std::cout << fragmentWeightMatrix.at(i * spectrumPeaks + j) << " ";
			}
			else {
				std::cout << (int)fragmentWeightMatrix.at(i * spectrumPeaks + j) << " ";
			}
			if (fragmentWeightMatrix[i * spectrumPeaks + j] > max_num_of_peaks) {
				max_num_of_peaks = fragmentWeightMatrix[i * spectrumPeaks + j];
			} 
			optimal_fragment_position = std::pair<int, int>(i, j);
			optimal_fragment_positions.push_back(optimal_fragment_position);
			fragment_to_num_of_peaks[optimal_fragment_position] = fragmentWeightMatrix[i * spectrumPeaks + j];
		}
	}
	std::vector<std::pair<double, std::pair<int, int> > > actual_optimal_fragment_positions;
	std::cout << "The maximum number of peaks in a potential peptide from the spectrum is: " + std::to_string(max_num_of_peaks);
	std::vector<std::pair<int, int> > optimal_fragment_positions_without_peaks;
	transform(actual_optimal_fragment_positions.begin(), actual_optimal_fragment_positions.end(), std::back_inserter(optimal_fragment_positions_without_peaks), [&](std::pair<double, std::pair<int, int> > actual_optimal_fragment_position) { return actual_optimal_fragment_position.second; });
	std::vector<std::vector<std::string> > best_fragments;
	for (int i = 0; i < optimal_fragment_positions_without_peaks.size(); i++) {
		best_fragments.push_back(FragmentWeightMatrix::find_sequence_from_fragment_position(fragmentWeightMatrix, edges, edges_without_amino_acids, actual_optimal_fragment_positions.at(i).first, optimal_fragment_positions_without_peaks.at(i), spectrumPeaks, prefix_suffix_connectors.at(i)));
	}
	for (int i = 0; i < optimal_fragment_positions_without_peaks.size(); i++) {
		std::cout << "The optimal peptide is with prefix-suffix connector " << prefix_suffix_connectors.at(i) << " is " << best_fragments.at(i).at(0) << "." << std::endl;
	}
}

// Initialization methods

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

// pybind11-related matters

PYBIND11_MODULE(denovo_sequencing, handle) {
	handle.doc() = "This is the DeNovo sequencing algorithm.";
	handle.def("main", 
    []() {
        pybind11::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
        main();
    });

	py::class_<AminoAcids>(handle, "AminoAcids").def(py::init<>());
}
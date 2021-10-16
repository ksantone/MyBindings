#include "SpectrumGraph.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "AminoAcids.h"
#include <map>


float round_var(float var)
{
	// 37.66666 * 100 =3766.66
	// 3766.66 + .5 =3767.16    for rounding off value
	// then type cast to int so value is 3767
	// then divided by 100 so the value converted into 37.67
	float value = (int)(var * 10 + .5);
	return (float)value / 10;
}


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
				for (it = AminoAcids::amino_acids.begin(); it != AminoAcids::amino_acids.end(); it++) {
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

/*SpectrumGraph construct_new_spectrum_graph(std::vector<int> desirable_peaks, std::vector<std::pair<int, int>> desirable_edges, SpectrumGraph spectrumGraph) {
	std::vector<Peak> newSpectralPeaks;
	std::map<int, int> map_to_new_indices;
	for (int i = 0; i < desirable_peaks.size(); i++) {
		map_to_new_indices.insert(std::pair<int, int>(desirable_peaks.at(i), newSpectralPeaks.size()));
		newSpectralPeaks.push_back(spectrumGraph.spectralPeaks.at(desirable_peaks.at(i)));
	}
	spectrumGraph.spectralPeaks = newSpectralPeaks;
	std::map<std::pair<int, int>, std::string> newSpectralEdges;
	for (int i = 0; i < desirable_edges.size(); i++) {
		newSpectralEdges[std::pair<int, int>(map_to_new_indices.at(std::get<0>(desirable_edges.at(i))), map_to_new_indices.at(std::get<1>(desirable_edges.at(i))))] = "";
	}
	spectrumGraph.spectralEdges = newSpectralEdges;
	return spectrumGraph;
}*/

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

/*std::vector<int> SpectrumGraph::compute_descendants_of_peak(std::vector<std::pair<int, int>> edges_for_current_peak, int i) {
	std::vector<int> peaks_to_return;
	peaks_to_return.push_back(std::get<1>(edges_for_current_peak.at(0)));
	if (i == 0) {
		peaks_to_return.push_back(std::get<1>(edges_for_current_peak.at(1)));
	}
	return peaks_to_return;
}

std::map<std::pair<int, int>, std::string> SpectrumGraph::compute_edges_from_peak(std::vector<std::pair<int, int>> edges_for_current_peak, int i) {
	std::map<std::pair<int, int>, std::string> edges_to_return;
	edges_to_return[edges_for_current_peak.at(0)] = ;
	if (i == 0) {
		edges_to_return[edges_for_current_peak.at(1)];
	}
	return edges_to_return;
}

SpectrumGraph SpectrumGraph::filter_edges_and_peaks(SpectrumGraph spectrumGraph) {
	std::vector<int> desirable_peaks;
	std::vector<std::pair<int, int>> desirable_edges;
	std::vector<int> processed_peaks;
	std::vector<std::pair<int, int>> allSpectralEdges = spectrumGraph.spectralEdges;
	for (int i = 0; i < allSpectralEdges.size(); i++) {
		int current_peak = std::get<0>(allSpectralEdges.at(i));
		std::vector<std::pair<int, int>> edges_for_current_peak = allSpectralEdges;
		if ((current_peak == 0 && std::find(desirable_peaks.begin(), desirable_peaks.end(), current_peak) == desirable_peaks.end())) {
			desirable_peaks.push_back(0);
			edges_for_current_peak.erase(std::remove_if(edges_for_current_peak.begin(), edges_for_current_peak.end(), [&](std::pair<int, int> edge) {return std::get<0>(edge)!=current_peak;}), edges_for_current_peak.end());
			std::vector<int> resulting_peaks = compute_descendants_of_peak(edges_for_current_peak, i);
			std::vector<std::pair<int, int>> resulting_edges = compute_edges_from_peak(edges_for_current_peak, i);
			resulting_peaks.erase(std::remove_if(resulting_peaks.begin(), resulting_peaks.end(), [&](int peak) { return std::find(desirable_peaks.begin(), desirable_peaks.end(), peak) != desirable_peaks.end(); }), resulting_peaks.end());
			desirable_peaks.insert(desirable_peaks.end(), resulting_peaks.begin(), resulting_peaks.end());
			desirable_edges.insert(desirable_edges.end(), resulting_edges.begin(), resulting_edges.end());
		}
		else if ((current_peak != 0 && std::find(desirable_peaks.begin(), desirable_peaks.end(), current_peak) != desirable_peaks.end() &&
					std::find(processed_peaks.begin(), processed_peaks.end(), current_peak) == processed_peaks.end())) {
			edges_for_current_peak.erase(std::remove_if(edges_for_current_peak.begin(), edges_for_current_peak.end(), [&](std::pair<int, int> edge) {return std::get<0>(edge) != current_peak;}), edges_for_current_peak.end());
			std::vector<int> resulting_peaks = compute_descendants_of_peak(edges_for_current_peak, i);
			resulting_peaks.erase(std::remove_if(resulting_peaks.begin(), resulting_peaks.end(), [&](int peak) { return std::find(desirable_peaks.begin(), desirable_peaks.end(), peak) != desirable_peaks.end(); }), resulting_peaks.end());
			std::vector<std::pair<int, int>> resulting_edges = compute_edges_from_peak(edges_for_current_peak, i);
			desirable_peaks.insert(desirable_peaks.end(), resulting_peaks.begin(), resulting_peaks.end());
			desirable_edges.insert(desirable_edges.end(), resulting_edges.begin(), resulting_edges.end());
			processed_peaks.push_back(current_peak);
		}
	}
	SpectrumGraph newSpectrumGraph = construct_new_spectrum_graph(desirable_peaks, desirable_edges, spectrumGraph);
	return newSpectrumGraph;
}*/
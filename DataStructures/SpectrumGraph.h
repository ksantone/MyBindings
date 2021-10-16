#include <vector>
#include <string>
#include "Peak.h"
#include <list>
#include <map>

#ifndef SpectrumGraph_HEADER
#define SpectrumGraph_HEADER

class SpectrumGraph {
public:
	std::vector<Peak> spectralPeaks;
	std::map<std::pair<int, int>, std::string> spectralEdges;
	std::map<std::pair<int, int>, std::string> computeSpectralEdges(std::vector<Peak>);
	//SpectrumGraph filter_edges_and_peaks(SpectrumGraph);
	//std::vector<int> compute_descendants_of_peak(std::vector<std::pair<int, int>>, int);
	std::map<std::pair<int, int>, std::string> compute_edges_from_peak(std::vector<std::pair<int, int> >, int);
	SpectrumGraph(std::string);
	SpectrumGraph(std::vector<int>, std::vector<std::pair<std::string, std::pair<int, int> > >);
};

#endif
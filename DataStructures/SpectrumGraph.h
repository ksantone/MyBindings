#include <vector>
#include <string>
#include "Peak.h"
#include <list>
#include <map>

#ifndef SpectrumGraph_HEADER
#define SpectrumGraph_HEADER

class SpectrumGraph {
public:
	double precursorMass;
	double precursorCharge;
	std::vector<Peak> spectralPeaks;
	std::vector<Peak> unmodifiedFinalSpectralPeaks;
	std::map<std::pair<int, int>, std::vector<std::string>> spectralEdges;
	std::vector<std::pair<std::vector<std::string>, std::vector<std::string> > > computeMassDecomposition(std::pair<int, int>, std::vector<std::vector<Peak> >);
	std::map<std::pair<int, int>, std::vector<std::string>> computeSpectralEdges(std::vector<Peak>);
	//SpectrumGraph filter_edges_and_peaks(SpectrumGraph);
	//std::vector<int> compute_descendants_of_peak(std::vector<std::pair<int, int>>, int);
	std::map<std::pair<int, int>, std::string> spectralPrefixSuffixConnectors;
	std::map<std::pair<int, int>, std::string> computePrefixSuffixConnectors(std::vector<Peak>, std::map<std::pair<int, int>, std::vector<std::string>>, float);
	std::map<std::pair<int, int>, std::string> compute_edges_from_peak(std::vector<std::pair<int, int> >, int);
	SpectrumGraph(std::string, double, double);
	SpectrumGraph(std::vector<int>, std::vector<std::pair<std::string, std::pair<int, int> > >);
};

#endif
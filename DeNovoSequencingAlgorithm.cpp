#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>
#include <queue>
#include <set>
#include <chrono>
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
        {"G", 57.02146}, {"H", 137.05891}, {"(I/L)", 113.08406}, {"K", 128.09496}, {"M", 131.04049}, {"F", 147.06841}, {"P", 97.05276},
        {"S", 87.03203}, {"T", 101.04768}, {"W", 186.07931}, {"Y", 163.06333}, {"V", 99.06841} };
std::mutex mtx;
std::condition_variable current_cv;
std::condition_variable threads_completing_cv;
std::vector<bool> triangle_complete;
std::mutex m1, m2, m3;
int threads_executing = 0;
int rows_completed = 0;

std::ofstream MyFile("/usr/src/app/output_data.txt");
std::ofstream MyMZs("/usr/src/app/output_mz.txt");
std::ofstream MyPrefixSuffix("/usr/src/app/output_prefix_suffix.txt");
std::ofstream MyTesting("/usr/src/app/testing.txt");
std::ofstream MyBIons("/usr/src/app/b_ions.txt");
std::ofstream MyPeaks("/usr/src/app/peaks.txt");
std::ofstream MyPreviousSequence("/usr/src/app/my_previous_sequence.txt");
std::ofstream MyFinalSequence("/usr/src/app/my_final_sequence.txt");

FILE *fp;

class AminoAcids {
public:
    static std::map <std::string, float> amino_acids;
    std::map<std::string, float> getAminoAcids();
    float getAminoAcidMass(std::string amino_acid);
    AminoAcids();
};

AminoAcids::AminoAcids() {

}

struct less_than_mz {
	inline bool operator()(const Peak& peak_1, const Peak& peak_2)
	{
		return (peak_1.mass < peak_2.mass);
	}
};

float round_var(float var)
{
	float value = (int)(var * 10 + .5);
	return (float)value / 10;
}

// SpectrumGraph-related methods

std::vector<Peak> glue_reversal(std::vector<Peak> tempSpectralPeaks, int charge, float max_mz) {
	tempSpectralPeaks.push_back(Peak(max_mz, charge, tempSpectralPeaks.at(tempSpectralPeaks.size()-1).intensity, tempSpectralPeaks.at(0).rt));
	//std::cout << "After gluing max mz size is: " << tempSpectralPeaks.size() << "." << std::endl;
	std::vector<Peak> finalSpectralPeaks = tempSpectralPeaks;
	for (int i = 1; i < tempSpectralPeaks.size()-1; i++) {
		//std::cout << "On peak: " << i << "." << std::endl;
		Peak tempPeak = tempSpectralPeaks.at(i);
		finalSpectralPeaks.push_back(Peak(max_mz-tempPeak.mass, charge, tempPeak.intensity, tempPeak.rt));
	}
	//std::cout << "After for loop size is: " << finalSpectralPeaks.size() << "." << std::endl;
	return finalSpectralPeaks;
}

std::map<std::pair<int, int>, std::string> SpectrumGraph::computePrefixSuffixConnectors(std::vector<Peak> spectralPeaks, std::map<std::pair<int, int>, std::string> spectralEdges, float maxpoint) {
	MyPrefixSuffix << "Maxpoint is: " << maxpoint << std::endl;
	//std::cout << "Peak at 45th index: " << spectralPeaks.at(44).mz << std::endl;
	std::map<std::pair<int, int>, std::string> prefix_suffix_connectors;
	int count = 1;
	/*while (spectralPeaks.at(count).mz < midpoint-187) {
		count++;
	}
	count++;*/
	std::map<std::string, float>::iterator aminoIt;
	bool notThereYet = true;
	MyPrefixSuffix << "Number of peaks is: " << spectralPeaks.size() << "." << std::endl;
	while (count < spectralPeaks.size()) {
		std::vector<int> indicesOfHigherPeaks;
		std::vector<std::string> aminoAcidsVector;
		// Loop through amino acids and check if there is a corresponding point within a threshold of this plus that amino acid mass
		for (aminoIt = amino_acids.begin(); aminoIt != amino_acids.end(); aminoIt++) {
			std::vector<Peak>::iterator it = std::find_if (spectralPeaks.begin(), spectralPeaks.end(), [&](Peak peak){return std::abs(maxpoint-(peak.mass+spectralPeaks.at(count).mass+aminoIt->second))<1;});
			if (it!=spectralPeaks.end()) {
				int index = it - spectralPeaks.begin();

				indicesOfHigherPeaks.push_back(index);
				aminoAcidsVector.push_back(aminoIt->first);
			}
		}
		for (int i = 0; i < indicesOfHigherPeaks.size(); i++) {
			std::vector<Peak>::iterator it = std::find_if (spectralPeaks.begin(), spectralPeaks.end(), [&](Peak peak){return std::abs((peak.mass+spectralPeaks.at(indicesOfHigherPeaks.at(i)).mass)-maxpoint)<2;});
			// TODO: Count all modified indices not just first
			MyPrefixSuffix << "Prefix-suffix of: " << aminoAcidsVector.at(i) << " being added with indicies (" << count << ", " << indicesOfHigherPeaks.at(i) << ") and masses (" << spectralPeaks.at(count).mass << ", " << spectralPeaks.at(indicesOfHigherPeaks.at(i)).mass << ")." << std::endl;
			//if (count==86 && indicesOfHigherPeaks.at(i)==123) {
				//std::cout << "The mass is: " << it->mass << "." << std::endl;
				//std::cout << "In here..." << std::endl;
				prefix_suffix_connectors.insert(std::pair<std::pair<int, int>, std::string>(std::pair<int, int>(count, indicesOfHigherPeaks.at(i)), aminoAcidsVector.at(i)));
			//}
		}
		count++;
	}
	std::cout << "Prefix-suffix connectors size is: " << prefix_suffix_connectors.size() << "." << std::endl;
	std::cout << "Max peak: " << maxpoint << "." << std::endl;
	
	return prefix_suffix_connectors;
}

std::map<std::pair<int, int>, std::string> SpectrumGraph::computeSpectralEdges(std::vector<Peak> spectralPeaks) {
	std::map<std::pair<int, int>, std::string> output_edges;
	for (int i = 0; i < spectralPeaks.size(); i++) {
		// Check if spectralPeaks differences are in dictionary (maybe make dictionary round down to near 1/10^n then check if either equal
		// or less than actual value by 1/10^n
		for (int j = i + 1; j < spectralPeaks.size(); j++) {
			float peak_diff = spectralPeaks.at(j).mass - spectralPeaks.at(i).mass;
			if (peak_diff >= 187) {
				break;
			}
			else if (peak_diff <= 56) {
				continue;
			}
			else {
				std::map<std::string, float>::iterator it;
				for (it = amino_acids.begin(); it != amino_acids.end(); it++) {
					if (std::abs((it->second)-peak_diff)<1.0) {
						MyMZs << "(" << i << ", " << j << ")" << std::endl;
						MyMZs << "(" << spectralPeaks.at(i).mass << ", " << spectralPeaks.at(j).mass << ")" << std::endl; 
						MyMZs << "Amino acid is: " << it->first << std::endl << std::endl; 
						//output_edges.emplace_back(std::pair<Peak, Peak>(spectralPeaks.at(i), spectralPeaks.at(j)));
						output_edges[std::pair <int, int>(i, j)] = it->first;
						// Sometimes there are multiple amino acids that match (unless we reduce to say 0.4, but for now let's just
						// take the first one)
					}
				}
			}
		}
	}
	MyMZs << "End of current charge..." << std::endl;
	return output_edges;
}

SpectrumGraph::SpectrumGraph(std::string spectralData, double precursor_mz, double precursor_charge) {
	precursorMass = precursor_mz * precursor_charge;
	precursorCharge = precursor_charge;
	std::cout << "Spectral data: " << spectralData << "." << std::endl;
	std::istringstream iss(spectralData);
	std::string del = ",";
	std::string peak_info;
	std::stringstream ss;
	float mass, rt, intensity;
	iss >> rt;
	std::vector<Peak> baseSpectralPeaks;
	baseSpectralPeaks.push_back(Peak(0, 1, 1, rt));
	//std::cout << "About to enter while loop..." << std::endl;
	while (iss >> peak_info) {
		int start = 0;
		int end = peak_info.find(del);
		ss << peak_info.substr(start, end - start) << std::endl;
		mass = std::stof(ss.str());
		start = end + del.size();
		end = peak_info.size();
		ss.str("");
		ss.clear();
		ss << peak_info.substr(start, end - start) << std::endl;
		intensity = std::stof(ss.str());
		baseSpectralPeaks.push_back(Peak(mass, 1, intensity, rt));
		ss.str("");
		ss.clear();
	}
	std::vector<Peak> unmodified_spectralPeaks;
	unmodified_spectralPeaks.insert(unmodified_spectralPeaks.end(), baseSpectralPeaks.begin(), baseSpectralPeaks.end());
	baseSpectralPeaks = glue_reversal(baseSpectralPeaks, 1, precursorMass);
	std::sort(unmodified_spectralPeaks.begin(), unmodified_spectralPeaks.end(), less_than_mz());
	std::sort(baseSpectralPeaks.begin(), baseSpectralPeaks.end(), less_than_mz());
	MyPeaks << "Base spectral peaks..." << std::endl;
	for (int i = 0; i < baseSpectralPeaks.size(); i++) {
		MyPeaks << "Mass: " << baseSpectralPeaks.at(i).mass << ", Charge: " << baseSpectralPeaks.at(i).charge << "." << std::endl;
	}
	for (int i = 0; i < unmodified_spectralPeaks.size(); i++) {
		MyPeaks << "Mass (unmodified): " << unmodified_spectralPeaks.at(i).mass << ", Charge: " << unmodified_spectralPeaks.at(i).charge << "." << std::endl;
	}
	std::vector<Peak>::iterator specIt;
	float maxPeak = (*(baseSpectralPeaks.end()-1)).mass;
	int count = 0;
	std::cout << "Max peak is: " << maxPeak/2 << std::endl;
	/*for (specIt = baseSpectralPeaks.begin(); (*specIt).mass < maxPeak/2; specIt++) {
		std::cout << "Current peak is: " << (*specIt).mass << "." << std::endl;
		count++;
	}
	std::cout << "Count is: " << count << std::endl;*/
	//std::vector<Peak>::iterator middlePeak = baseSpectralPeaks.begin() + count;
	/*MyPeaks << "Middle peak is: " << (*middlePeak).mass << std::endl;
	std::vector<Peak> newBaseSpectralPeaks = std::vector<Peak>(baseSpectralPeaks.begin(), middlePeak);*/
	SpectrumGraph::unmodifiedFinalSpectralPeaks.insert(SpectrumGraph::unmodifiedFinalSpectralPeaks.end(), unmodified_spectralPeaks.begin(), unmodified_spectralPeaks.end());
	SpectrumGraph::spectralPeaks.insert(SpectrumGraph::spectralPeaks.end(), baseSpectralPeaks.begin(), baseSpectralPeaks.end());
	for (int i = 2; i <= precursorCharge; i++) {
		std::vector<Peak> tempSpectralPeaks;
		std::vector<Peak> tempUnmodified_SpectralPeaks;
		int count = 1;
		while (baseSpectralPeaks.at(count).mass<precursorMass/i) {
			//std::cout << "Current peak has MZ: " << baseSpectralPeaks.at(count).mz << std::endl;
			Peak tempPeak = Peak(baseSpectralPeaks.at(count).mass*i, i, baseSpectralPeaks.at(count).intensity, baseSpectralPeaks.at(count).rt);
			tempSpectralPeaks.push_back(tempPeak);
			count += 1;
		}
		count = 1;
		while (unmodified_spectralPeaks.at(count).mass<precursorMass/i) {
			//std::cout << "Current peak has MZ: " << baseSpectralPeaks.at(count).mz << std::endl;
			Peak tempUnmodifiedPeak = Peak(unmodified_spectralPeaks.at(count).mass*i, i, unmodified_spectralPeaks.at(count).intensity, unmodified_spectralPeaks.at(count).rt);
			tempUnmodified_SpectralPeaks.push_back(tempUnmodifiedPeak);
			count += 1;
		}
		std::sort(tempUnmodified_SpectralPeaks.begin(), tempUnmodified_SpectralPeaks.end(), less_than_mz());
		tempSpectralPeaks = glue_reversal(tempSpectralPeaks, i, precursorMass);
		std::sort(tempSpectralPeaks.begin(), tempSpectralPeaks.end(), less_than_mz());
		MyPeaks << "Spectral peaks with charge: " << i << "." << std::endl;
		for (int i = 0; i < tempSpectralPeaks.size(); i++) {
			MyPeaks << "Mass: " << tempSpectralPeaks.at(i).mass << ", Charge: " << tempSpectralPeaks.at(i).charge << "." << std::endl;
		}
		for (int i = 0; i < unmodified_spectralPeaks.size(); i++) {
			MyPeaks << "Mass (unmodified): " << unmodified_spectralPeaks.at(i).mass << ", Charge: " << unmodified_spectralPeaks.at(i).charge << "." << std::endl;
		}
		std::cout << "AAA" << std::endl;
		/*std::vector<Peak>::iterator specIt;
		float maxPeak = (*(tempSpectralPeaks.end()-1)).mass;
		std::cout << "BBB" << std::endl;
		for (specIt = tempSpectralPeaks.begin(); (*specIt).mass < maxPeak/2; specIt++) {
			std::cout << "Current mass is: " << (*specIt).mass << " and count is: " << count+1 << std::endl;
			count++;
		}
		std::vector<Peak>::iterator middlePeak = tempSpectralPeaks.begin() + count;
		std::vector<Peak> newTempSpectralPeaks = std::vector<Peak>(tempSpectralPeaks.begin(), middlePeak);
		std::cout << "Number of peaks in charge " << i << ": " << newTempSpectralPeaks.size() << "." << std::endl;*/
		SpectrumGraph::unmodifiedFinalSpectralPeaks.insert(SpectrumGraph::unmodifiedFinalSpectralPeaks.end(), tempUnmodified_SpectralPeaks.begin(), tempUnmodified_SpectralPeaks.end());
		SpectrumGraph::spectralPeaks.insert(SpectrumGraph::spectralPeaks.end(), tempSpectralPeaks.begin(), tempSpectralPeaks.end());
		std::sort(unmodifiedFinalSpectralPeaks.begin(), unmodifiedFinalSpectralPeaks.end(), less_than_mz());
		std::sort(SpectrumGraph::spectralPeaks.begin(), SpectrumGraph::spectralPeaks.end(), less_than_mz());
	}
	MyPeaks << "Final peaks..." << std::endl;
	for (int i = 0; i < SpectrumGraph::spectralPeaks.size(); i++) {
		MyPeaks << "Mass: " << SpectrumGraph::spectralPeaks.at(i).mass << ", Charge: " << SpectrumGraph::spectralPeaks.at(i).charge << "." << std::endl;
	}
	for (int i = 0; i < SpectrumGraph::unmodifiedFinalSpectralPeaks.size(); i++) {
		MyPeaks << "Mass (unmodified): " << SpectrumGraph::unmodifiedFinalSpectralPeaks.at(i).mass << ", Charge: " << SpectrumGraph::unmodifiedFinalSpectralPeaks.at(i).charge << "." << std::endl;
	}
	std::cout << "Before computeSpectralEdges..." << std::endl;
	SpectrumGraph::spectralEdges = computeSpectralEdges(SpectrumGraph::spectralPeaks);
	std::cout << "After computeSpectralEdges..." << std::endl;
	std::cout << "Before computePrefixSuffixConnectors..." << std::endl;
	spectralPrefixSuffixConnectors = computePrefixSuffixConnectors(SpectrumGraph::unmodifiedFinalSpectralPeaks/*SpectrumGraph::spectralPeaks*/, SpectrumGraph::spectralEdges, maxPeak);
	std::cout << "After computePrefixSuffixConnectors..." << std::endl;
	MyMZs.close();
	std::cout << "Number of spectral peaks vectors is: " << SpectrumGraph::spectralPeaks.size() << "." << std::endl;
	std::cout << "Number of spectral edges vectors is: " << SpectrumGraph::spectralEdges.size() << "." << std::endl;
}

Peak::Peak(float mass, int charge, float intensity, float rt)
{
	this->mass = mass;
	this->charge = charge;
	this->intensity = intensity;
	this->rt = rt;
}

// FragmentWeightMatrix-related methods

double FragmentWeightMatrix::findMaximumSubfragment(std::vector<double> fragmentWeightMatrix, int i, int j, std::map<std::pair<int, int>, std::string> edges, int spectrumPeaks) {
	double currentMaxfragment = -1;
	std::vector<std::pair<int, int> > edges_without_amino_acids;
	std::transform(edges.begin(), edges.end(), std::back_inserter(edges_without_amino_acids), [&](std::pair<std::pair<int, int>, std::string> edge) { return edge.first; });
	std::pair<int, int> final_pair;
	std::string amino_acid;
	for (int k = j - 1; k >= 0; k--) {
		/*if(i==0 && k==0) {
			std::cout << "At (0, 0) and value of j is: " << j << "." << std::endl;
			for (int i = 0; i < edges_without_amino_acids.size(); i++) {
				std::cout << "(" << edges_without_amino_acids.at(i).first << ", " <<edges_without_amino_acids.at(i).second << ")." << std::endl;
			}
		}*/
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
		//std::cout << "Computed fragment at position: (" << i << ", " << j << ") to have " << currentMaxfragment << " peaks." << std::endl;
		mtx.unlock();
	}
	return currentMaxfragment;
}

std::vector<std::vector<std::pair<int, int> > > FragmentWeightMatrix::glue_edge_to_position_vector(std::pair<int, int> current_edge, std::vector<std::vector<std::pair<int, int> > > edge_vectors, bool left) {
	std::vector<std::vector<std::pair<int, int> > > new_edge_vectors;
	for (int i = 0; i < edge_vectors.size(); i++) {
		std::vector<std::pair<int, int> > temp_vector = edge_vectors.at(i);
		temp_vector.push_back(current_edge);
		if (left) {
			if (edge_vectors.at(i).size()>0 && current_edge.second==edge_vectors.at(i).at(0).first) {
				std::rotate(temp_vector.rbegin(), temp_vector.rbegin()+1, temp_vector.rend());
				new_edge_vectors.push_back(temp_vector);
			}
		}
		else {
			if (edge_vectors.at(i).size()>0 && current_edge.first==edge_vectors.at(i).at(edge_vectors.at(i).size()-1).second) {
				new_edge_vectors.push_back(temp_vector);
			}
		}
	}
	return new_edge_vectors;
}

std::vector<std::string> FragmentWeightMatrix::glue_edge_to_all_sequences(std::string amino_acid, std::pair<int, int> current_edge, std::vector<std::string> subpeptides, std::vector<std::vector<std::pair<int, int> > > edge_vectors, bool left) {
	std::vector<std::string> new_subpeptides;
	for (int i = 0; i < subpeptides.size(); i++) {
		if (left) {
			if (edge_vectors.at(i).size()>0 && current_edge.second==edge_vectors.at(i).at(0).first) {
				new_subpeptides.push_back(amino_acid + subpeptides.at(i));
			}
		}
		else {
			if (edge_vectors.at(i).size()>0 && current_edge.first==edge_vectors.at(i).at(edge_vectors.at(i).size()-1).second) {
				new_subpeptides.push_back(subpeptides.at(i) + amino_acid);
			}
		}
	}
	return new_subpeptides;
}

std::vector<std::vector<std::pair<std::string, int> > > computeMassDecomposition(float ion_mass, int amino_index, std::vector<std::pair<std::string, float> > amino_acids_vector) {
	std::vector<std::vector<std::pair<std::string, int> > > mass_decompositions;
	if (amino_index>0) {
		//MyBIons << amino_acids_vector.at(amino_index-1).first << " and ion mass is " << ion_mass << "." << std::endl;
		/*if(amino_acids_vector.at(amino_index-1).first=="E") {
			if(std::abs(ion_mass-454.24)<0.1) {
				MyTesting << "E" << std::endl;
			}
		}
		if(amino_acids_vector.at(amino_index-1).first=="H") {
			if(std::abs(ion_mass-317.181)<0.1) {
				MyTesting << "EH" << std::endl;
			}
		}
		if(amino_acids_vector.at(amino_index-1).first=="K") {
			if(std::abs(ion_mass-189.086)<0.1) {
				MyTesting << "EHK" << std::endl;
			}
		}
		if(amino_acids_vector.at(amino_index-1).first=="S") {
			if(std::abs(ion_mass-102.054)<0.1) {
				MyTesting << "EHKS" << std::endl;
			}
		}
		if(amino_acids_vector.at(amino_index-1).first=="T") {
			if(std::abs(ion_mass)<2.0) {
				MyTesting << "EHKST" << std::endl;
			}
		}*/
	}
	if (amino_index!=19) {
		for (int i = 0; i < std::min(3, 1+(int)(ion_mass/amino_acids_vector.at(amino_index).second)); i++) {
			if (i==0 && amino_acids_vector.at(amino_index).first=="T" && std::abs(ion_mass-102.054)<0.1) {
				MyTesting << "Ion mass is: " << ion_mass << "." << std::endl;
				MyTesting << "Getting somewhere (EHKST): " << ion_mass-amino_acids_vector.at(amino_index).second << "." << std::endl;
			}
			if (i==1 && amino_acids_vector.at(amino_index).first=="T" && std::abs(ion_mass-102.054)<0.1) {
				MyTesting << "Getting somewhere (EHKST): " << ion_mass-i*amino_acids_vector.at(amino_index).second << "." << std::endl;
			}
			std::vector<std::vector<std::pair<std::string, int> > > new_mass_decompositions;
			if (ion_mass-i*amino_acids_vector.at(amino_index).second>2.0) {
				new_mass_decompositions = computeMassDecomposition(ion_mass-i*amino_acids_vector.at(amino_index).second, amino_index+1, amino_acids_vector);
				for (int j = 0; j < new_mass_decompositions.size(); j++) {
					new_mass_decompositions.at(j).push_back(std::pair<std::string, int>(amino_acids_vector.at(amino_index).first, i));
				}
			} else {
				std::vector<std::pair<std::string, int> > new_mass_decomposition;
				new_mass_decomposition.push_back(std::pair<std::string, int>(amino_acids_vector.at(amino_index).first, i));
				new_mass_decompositions.push_back(new_mass_decomposition);
			}
			mass_decompositions.insert(mass_decompositions.end(), new_mass_decompositions.begin(), new_mass_decompositions.end());
		}
	}
	return mass_decompositions;
}

/*std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > FragmentWeightMatrix::find_sequence_from_fragment_position(std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > final_sequences, std::vector<double> fragmentWeightMatrix, std::map<std::pair<int, int>, std::string> edges, std::vector<std::pair<int, int> > edges_without_amino_acids, int spectrumPeaks, std::vector<std::pair<int, int> > current_edges, std::string prefix_suffix_connector) {
	if (final_sequences.first.size()==0) {
		final_sequences.first.push_back(prefix_suffix_connector);
		final_sequences.second.push_back(current_edges);
		std::cout << "Prefix-suffix connector is: " << prefix_suffix_connector << "." << std::endl;
		//MyPrefixSuffix << "Prefix-suffix connector is: " << prefix_suffix_connector << "." << std::endl;
		//MyPrefixSuffix << "Corresponding to edge: (" << current_edges.at(0).first << ", " << current_edges.at(0).second << ")." << std::endl;
	}
	//double max_prev_frag = -std::numeric_limits<double>::infinity();
	//std::cout << "Current edge: (" << current_edges.at(0).first << ", " << current_edges.at(0).second << ")." << std::endl;
	std::queue<int> left_edges, right_edges;
	left_edges.push(current_edges.at(0).first);
	right_edges.push(current_edges.at(0).second);
	int i, j;
	std::pair<int, int> prev_frag_pair;
	bool found;
	std::vector<std::string> vector_1, vector_2;
	std::vector<std::vector<std::pair<int, int> > > position_vector_1, position_vector_2;
	std::set<int> left_used, right_used;
	while(true) {
		found = false;
		if (!left_edges.empty()) {
			i = left_edges.front();
			left_edges.pop();
		} else {
			i = 0;
		}
		left_used.insert(i);
		//std::cout << "About to enter first loop with i = " << i << "." << std::endl;
		// Use 2 queue to push elements in (when if statement is hit) then set i and j to the first values in the queue until the queue is empty
		for (int k = i - 1; k >= 0; k--) {
			if (std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, i))) {
				//max_prev_frag = fragmentWeightMatrix.at(k * spectrumPeaks + j);
				prev_frag_pair = std::pair<int, int>(k, i);
				found = true;
				vector_1 = FragmentWeightMatrix::glue_edge_to_all_sequences(edges[prev_frag_pair], prev_frag_pair, final_sequences.first, final_sequences.second, true);
				position_vector_1 = FragmentWeightMatrix::glue_edge_to_position_vector(prev_frag_pair, final_sequences.second, true);
				final_sequences.first.insert(final_sequences.first.end(), vector_1.begin(), vector_1.end());
				final_sequences.second.insert(final_sequences.second.end(), position_vector_1.begin(), position_vector_1.end());
				if (k!=0 && left_used.find(k)==left_used.end()) {
					left_used.insert(k);
					left_edges.push(k);
				}
			}
		}
		if (!right_edges.empty()) {
			j = right_edges.front();
			right_edges.pop();
		} else {
			j = 0;
		}
		right_used.insert(j);
		//std::cout << "About to enter first loop with j = " << j << "." << std::endl;
		for (int k = j - 1; k >= 0; k--) {
			if (std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, j))) {
				//max_prev_frag = fragmentWeightMatrix.at(i * spectrumPeaks + k);
				prev_frag_pair = std::pair<int, int>(k, j);
				found = true;
				//std::cout << "Size of first: " << final_sequences.first.size() << "." << std::endl;
				//std::cout << "Size of second: " << final_sequences.second.size() << "." << std::endl;
				vector_2 = FragmentWeightMatrix::glue_edge_to_all_sequences(edges[prev_frag_pair], prev_frag_pair, final_sequences.first, final_sequences.second, false);
				position_vector_2 = FragmentWeightMatrix::glue_edge_to_position_vector(prev_frag_pair, final_sequences.second, false);
				final_sequences.first.insert(final_sequences.first.end(), vector_2.begin(), vector_2.end());
				final_sequences.second.insert(final_sequences.second.end(), position_vector_2.begin(), position_vector_2.end());
				if (k!=0 && right_used.find(k)==right_used.end()) {
					right_used.insert(k);
					right_edges.push(k);
				}
			}
		}
		if(i==0 && j==0) {
			break;
		}
		if (i==0 && j<10) {
			int count = left_edges.size();
			while(!left_edges.empty()) {
				std::cout << left_edges.front() << std::endl;
				left_edges.pop();
			}
			break;
		}
		// Change exit condition to both queues being empty
	}
	return final_sequences;
}*/

std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > FragmentWeightMatrix::find_b_sequence_from_fragment_position(std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > final_sequences, std::vector<double> fragmentWeightMatrix, std::map<std::pair<int, int>, std::string> edges, std::vector<std::pair<int, int> > edges_without_amino_acids, int spectrumPeaks, std::vector<std::pair<int, int> > current_edges, std::string prefix_suffix_connector, std::vector<std::vector<std::pair<std::string, int> > > mass_decompositions, std::vector<Peak> unmodifiedFinalSpectralPeaks) {
	MyFinalSequence << "First peak has mass: " << unmodifiedFinalSpectralPeaks.at(0).mass << "." << std::endl;
	if (final_sequences.first.size()==0) {
		std::cout << "Inside..." << std::endl;
		//final_sequences.first.push_back(prefix_suffix_connector);
		final_sequences.second.push_back(current_edges);
		std::cout << "Prefix-suffix connector is: " << prefix_suffix_connector << "." << std::endl;
	}
	std::cout << "After if..." << std::endl;
	std::set<int> left_edges, right_edges;
	// Should be first, but for some reason need second...?
	MyFinalSequence << "Maximum mass of b-ions is: " << unmodifiedFinalSpectralPeaks.at(current_edges.at(0).first).mass << "." << std::endl;
	left_edges.insert(current_edges.at(0).first);
	right_edges.insert(0);
	std::cout << "After pushes..." << std::endl;
	int i;
	std::vector<std::string> vector;
	std::vector<std::vector<std::pair<int, int> > > position_vector;
	std::set<int> used;
	bool found;
	int min_peak;
	std::string current_amino;
	std::vector<std::pair<std::string, int> > fixed_mass_decomposition = mass_decompositions.at(0); // Will unfix later
	std::vector<std::string> fixed_mass_decomp;
	for (int i = 0; i < fixed_mass_decomposition.size(); i++) {
		if (fixed_mass_decomposition.at(i).second>0) {
			fixed_mass_decomp.push_back(fixed_mass_decomposition.at(i).first);
		}
	}
	MyFinalSequence << "Size of mass decomposition is: " << fixed_mass_decomp.size() << "." << std::endl;
	MyFinalSequence << "Spectral peaks size is: " << unmodifiedFinalSpectralPeaks.size() << "." << std::endl;
	std::vector<std::pair<std::string, int> > fixed_mass_decomposition_pairs;
	std::vector<std::pair<float, std::string> > b_ion_prefix_tree;
	std::queue<std::string> unsearched_sequences;
	std::map<std::string, int> sequence_to_error;
	std::map<std::string, int> sequence_to_peak;
	b_ion_prefix_tree.push_back(std::pair<float, std::string>(0, ""));
	unsearched_sequences.push("");
	sequence_to_error[""] = 0;
	sequence_to_peak[""] = 0;
	std::cout << "Just before single or at least pair while loop..." << std::endl;
	MyFinalSequence << "Searching for prefixes now..." << std::endl;
	std::string previous_seq;
	while (true) {
		while (!unsearched_sequences.empty()) {
			previous_seq = unsearched_sequences.front();
			unsearched_sequences.pop();
			if (sequence_to_error[previous_seq]<=1) {
				break;
			}
		}
		MyFinalSequence << "Previous sequence is: " << previous_seq << "." << std::endl;
		float sequence_mass = 0;
		for (int k = 0; k < previous_seq.length(); k++) {
			sequence_mass += amino_acids[std::string(1, previous_seq.at(k))];
		}
		float temp_sequence_mass = sequence_mass;
		for (int j = 0; j < fixed_mass_decomp.size(); j++) {
			//MyFinalSequence << "On amino acid: " << fixed_mass_decomp.at(j) << "." << std::endl;
			current_amino = fixed_mass_decomp.at(j);
			sequence_mass = temp_sequence_mass + amino_acids[current_amino];
			bool found_first = false, sufficient_distance = false;
			int peak_index = sequence_to_peak[previous_seq];
			//MyFinalSequence << "Peak index is: " << peak_index << std::endl;
			while(peak_index<=current_edges.at(0).first) {
				if (std::abs(unmodifiedFinalSpectralPeaks.at(peak_index).mass-sequence_mass)<17) {
					if (std::abs(unmodifiedFinalSpectralPeaks.at(peak_index).mass-sequence_mass)<3) {
						sufficient_distance = true;
					}
					if (!std::count(b_ion_prefix_tree.begin(), b_ion_prefix_tree.end(), std::pair<float, std::string>(unmodifiedFinalSpectralPeaks.at(peak_index).mass, previous_seq+current_amino))) {
						sequence_to_peak[previous_seq+current_amino] = peak_index;
						unsearched_sequences.push(previous_seq+current_amino);
						b_ion_prefix_tree.push_back(std::pair<float, std::string>(unmodifiedFinalSpectralPeaks.at(peak_index).mass, previous_seq+current_amino));
						MyFinalSequence << "Mass: " << unmodifiedFinalSpectralPeaks.at(peak_index).mass << " partial sequence: " << previous_seq+current_amino << "." << std::endl;
						found_first = true;
					}
					peak_index += 1;
				} else if(found_first) {
					break;
				} else {
					peak_index += 1;
				}
			}
			MyFinalSequence << "Sufficient distance is: " << sufficient_distance << " on sequence " << previous_seq+current_amino << "." << std::endl;
			if (!sufficient_distance) {
				sequence_to_error[previous_seq+current_amino] = sequence_to_error[previous_seq]+1;
			}
		}
		if (unsearched_sequences.empty()) {
			break;
		}
	}
	for (int i = 0; i < b_ion_prefix_tree.size(); i++) {
		// Only include prefixes that go up to prefix-suffix connector and put into separate file
		//if (b_ion_prefix_tree.at(i).first==583.283) {
			MyFinalSequence << "Sequence at mass: " << b_ion_prefix_tree.at(i).first << " is: " << b_ion_prefix_tree.at(i).second << "." << std::endl;
		//}
	}
	//MyFinalSequence << "Final prefix sequence is: " << final_sequence_prefix << "." << std::endl;
	//MyFinalSequence << "Final suffix sequence is: " << final_sequence_suffix << "." << std::endl;
	//MyFinalSequence.close();
	return final_sequences;
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

FragmentWeightMatrix::FragmentWeightMatrix(SpectrumGraph spectrumGraph, int spectrumPeaks, std::map<std::pair<int, int>, std::string> edges) {
	std::cout << "In FragmentWeightMatrix constructor..." << std::endl;
	FragmentWeightMatrix::spectrumPeaks = spectrumPeaks;
	FragmentWeightMatrix::edges = edges;
	fp = fopen("/usr/src/app/progress.txt", "w+");
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
		//std::cout << "On thread number: " << i+1 << "." << std::endl;
		traverse_thread = std::thread([&]() {this->FragmentWeightMatrix::compute_row(i);});
		traverse_thread.join();
		//std::cout << "Thread number " << i+1 << " completed." << std::endl;
		fprintf(fp, "%d \n", i);
	}
	pybind11::gil_scoped_acquire acquire;
	double max_num_of_peaks = 0;
	std::pair<int, int> optimal_fragment_position;
	std::vector<std::pair<int, int> > edges_without_amino_acids;
	std::transform(edges.begin(), edges.end(), std::back_inserter(edges_without_amino_acids), [&](std::pair<std::pair<int, int>, std::string> edge) { return edge.first; });
	std::vector<std::string> prefix_suffix_connectors;
	std::vector<std::pair<int, int> > optimal_fragment_positions;
	std::map<std::pair<int, int>, double> fragment_to_num_of_peaks;
	for (int i = 0; i < spectrumPeaks; i++) {
		for (int j = 0; j < spectrumPeaks; j++) {
			if (fragmentWeightMatrix[i * spectrumPeaks + j] > max_num_of_peaks) {
				max_num_of_peaks = fragmentWeightMatrix[i * spectrumPeaks + j];
			}
			if (fragmentWeightMatrix[i * spectrumPeaks + j] > -1) {
				optimal_fragment_position = std::pair<int, int>(i, j);
				optimal_fragment_positions.push_back(optimal_fragment_position);
				fragment_to_num_of_peaks[optimal_fragment_position] = fragmentWeightMatrix[i * spectrumPeaks + j];
			}
		}
	}
	std::cout << "The maximum number of peaks in a potential peptide from the spectrum is: " + std::to_string(max_num_of_peaks);
	std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > > best_fragments;
	std::cout << "Size of current prefix-suffix connector iterator is: " << spectrumGraph.spectralPrefixSuffixConnectors.size() << "." << std::endl;
	std::vector<std::pair<int, int> > standing_edges;
	std::string prefix_suffix;
	std::vector<std::pair<float, std::vector<std::vector<std::pair<std::string, int> > > > > b_ion_to_mass_decompositions;
	bool first = true;
	for (std::map<std::pair<int, int>, std::string>::iterator it = spectrumGraph.spectralPrefixSuffixConnectors.begin(); it != spectrumGraph.spectralPrefixSuffixConnectors.end(); it++) {
		if (!first) {
			break;
		}
		std::cout << "Prefix-suffix connector is: " << it->second << "." << std::endl;
		std::cout << "Prefix-suffix connector edge is: (" << it->first.first << ", " << it->first.second << ")." << std::endl;
		std::cout << "Prefix-suffix connector edge masses are: (" << spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass << ", " << spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.second).mass << ")." << std::endl;
		std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > final_sequences;
		std::vector<std::pair<int, int>> current_edges;
		current_edges.push_back(it->first);
		standing_edges.push_back(it->first);
		prefix_suffix = it->second;
		//float ion_mass = spectrumGraph.spectralPeaks.at(amino_acids[it->second]);
		//std::cout << "Ion mass is: " << ion_mass << "." << std::endl;
		// Determine ion_mass from it->first and spectrumGraph.spectralPeaks (subtract mass of hydrogen=1.007825Da from b_ion_mass)
		std::vector<std::pair<std::string, float> > amino_acids_vector(amino_acids.begin(), amino_acids.end());
		std::cout << "Just before computation of mass decomposition." << std::endl;
		MyBIons << "Current mass is: " << spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass << " and truth value is: " << (std::abs(spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass-583.283)>=0.01) << "." << std::endl;
		if (std::abs(spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass-583.283)>=0.01) {
			continue;
		} else {
			MyBIons << "Inside of here..." << std::endl;
		}
		std::vector<std::vector<std::pair<std::string, int> > > mass_decompositions = computeMassDecomposition(spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass, 0, amino_acids_vector);
		MyBIons << "Mass decompositions computed for 583.283..." << std::endl;
		if (mass_decompositions.size()==0) {
			MyFinalSequence << "Empty mass decomposition." << std::endl;
			continue;
		}
		std::vector<std::vector<std::pair<std::string, int> > > temp_mass_decompositions;
		for (int i = 0; i < mass_decompositions.size(); i++) {
			bool K = false, T = false, E = false, S = false, H = false;
			for (int j = 0; j < mass_decompositions.at(i).size(); j++) {
				if (mass_decompositions.at(i).at(j).first=="K") {
					if (mass_decompositions.at(i).at(j).second>0) {
						K = true;
					} else {
						break;
					}
				} else if (mass_decompositions.at(i).at(j).first=="T") {
					if (mass_decompositions.at(i).at(j).second>0) {
						T = true;
					} else {
						break;
					}
				} else if (mass_decompositions.at(i).at(j).first=="E") {
					if (mass_decompositions.at(i).at(j).second>0) {
						E = true;
					} else {
						break;
					}
				} else if (mass_decompositions.at(i).at(j).first=="S") {
					if (mass_decompositions.at(i).at(j).second>0) {
						S = true;
					} else {
						break;
					}
				} else if (mass_decompositions.at(i).at(j).first=="H") {
					if (mass_decompositions.at(i).at(j).second>0) {
						H = true;
					} else {
						break;
					}
				} else {
					if (mass_decompositions.at(i).at(j).second>0) {
						K = false;
						break;
					}
				}
			}
			if (K==true && T==true && E==true && S==true && H==true) {
				temp_mass_decompositions.push_back(mass_decompositions.at(i));
			}
		}
		MyBIons << "Number of mass decomposition containing KTESH is: " << temp_mass_decompositions.size() << "." << std::endl;
		b_ion_to_mass_decompositions.push_back(std::pair<float, std::vector<std::vector<std::pair<std::string, int> > > >(spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass, mass_decompositions));
		MyBIons << "Number of mass decompositions for mass " << spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass << " is: " << mass_decompositions.size() << "." << std::endl;
		for (int i = 0; i < temp_mass_decompositions.size(); i++) {
				for(int j = 0; j < temp_mass_decompositions.at(i).size(); j++) {
					MyBIons << "(" << temp_mass_decompositions.at(i).at(j).first << ", " << temp_mass_decompositions.at(i).at(j).second << "), ";
				}
				MyBIons << std::endl;
		}
		MyBIons << "Now all mass decompositions total..." << std::endl;
		for (int i = 0; i < mass_decompositions.size(); i++) {
				for(int j = 0; j < mass_decompositions.at(i).size(); j++) {
					MyBIons << "(" << mass_decompositions.at(i).at(j).first << ", " << mass_decompositions.at(i).at(j).second << "), ";
				}
				MyBIons << std::endl;
		}
		//b_ion_to_mass_decompositions[spectrumGraph.spectralPeaks.at(it->first.first).mass] = mass_decompositions;
		MyFinalSequence << "At mass: " << spectrumGraph.unmodifiedFinalSpectralPeaks.at(it->first.first).mass << std::endl;
		std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > b_fragment = FragmentWeightMatrix::find_b_sequence_from_fragment_position(final_sequences, fragmentWeightMatrix, edges, edges_without_amino_acids, spectrumPeaks, current_edges, /*prefix_suffix*/"K", temp_mass_decompositions, spectrumGraph.unmodifiedFinalSpectralPeaks);
		// Make 2 separate find_b_sequence_from_fragment_position and find_y_sequence_from_fragment_position functions for finding peptide sequences
		//std::pair<std::vector<std::string>, std::vector<std::vector<std::pair<int, int> > > > new_fragment = FragmentWeightMatrix::find_sequence_from_fragment_position(final_sequences, fragmentWeightMatrix, edges, edges_without_amino_acids, spectrumPeaks, current_edges, it->second);
		//best_fragments.push_back(new_fragment);
		first = false;
	}
	MyBIons.close();
	MyFinalSequence.close();
	MyPreviousSequence.close();
	// Update class to have find_b_sequence_from_fragment_position and find_y_sequence_from_fragment_position
	//MyFinalSequence << "About to start computing final sequence..." << std::endl;
	std::cout << "Optimal fragments size: " << optimal_fragment_positions.size() << std::endl;
	// TODO don't take only first element in best fragments per element
	for (int j = 0; j < best_fragments.size(); j++) {
		for (int k = 0; k < best_fragments.at(j).first.size(); k++) {
			MyFile << best_fragments.at(j).first.at(k) << std::endl;
			for (int l = 0; l < best_fragments.at(j).second.at(k).size(); l++) {
				MyFile << "(" << best_fragments.at(j).second.at(k).at(l).first << ", " << best_fragments.at(j).second.at(k).at(l).second << ") ";
			}
			MyFile << std::endl;
		}
	}
	MyFile << "End of current charge..." << std::endl;
}

// Initialization methods

void DeNovoSequencingAlgorithm::run_denovo() {
    std::string spectraFileName = "/usr/src/app/new_spectrum_list.txt";
    std::cout << "File location is: " << spectraFileName << std::endl;
    std::ifstream spectraFile(spectraFileName);
    bool mzAndCharge = false, onPeaks = false;
    float precursor_mz, precursor_charge;
    std::list<SpectrumGraph> spectrumGraphs, filteredSpectrumGraphs;
    std::string spectralData;
    int count = 1;
    while (getline(spectraFile, spectralData)) {
    	std::cout << "On line " << count << std::endl;
    	std::cout << "Line is: " << spectralData << std::endl;
    	if (mzAndCharge) {
			int end = spectralData.find(",");
    		precursor_mz = std::stof(spectralData.substr(0, end));
    		precursor_charge = std::stof(spectralData.substr(end+1, spectralData.size()));
    		std::cout << "MZ: " << precursor_mz << " Charge: " << precursor_charge << std::endl;
    	} else if (onPeaks) {
    		spectrumGraphs.push_back(SpectrumGraph(spectralData, precursor_mz, precursor_charge));
    	}
        if (spectralData.find("Precursor")==0) {
        	mzAndCharge = true;
        	onPeaks = false;
        } else if (spectralData.find("RT")==0) {
        	std::cout << "About to call spectrum graph constructor..." << std::endl;
        	mzAndCharge = false;
        	onPeaks = true;
        } else {
        	mzAndCharge = false;
        	onPeaks = false;
        }
    }
    std::list<SpectrumGraph>::iterator it;
    std::list<FragmentWeightMatrix> fragmentWeightMatrices;
    std::cout << "Spectrum graphs size is: " << spectrumGraphs.size() << "." << std::endl;
    for (it = spectrumGraphs.begin(); it != spectrumGraphs.end(); it++) {
    	std::cout << "About to call fragment weight matrix constructor..." << std::endl;
        fragmentWeightMatrices.push_back(FragmentWeightMatrix(*it, it->spectralPeaks.size(), it->spectralEdges));
    }
    MyFile.close();
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
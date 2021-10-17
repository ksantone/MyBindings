#include <map>
#include "AminoAcids.h"

std::map <std::string, float> amino_acids = { {"A", 71.03711}, {"R", 156.10111}, {"N", 114.04293}, {"D", 115.02694}, {"C", 103.00919}, {"E", 129.04259}, {"Q", 128.05858 },
        {"G", 57.02146}, {"H", 137.05891}, {"I/L", 113.08406}, {"K", 128.09496}, {"M", 131.04049}, {"F", 147.06841}, {"P", 97.05276},
        {"S", 87.03203}, {"T", 101.04768}, {"W", 186.07931}, {"Y", 163.06333}, {"V", 99.06841} };

float getAminoAcidMass(std::string amino_acid) {
	return AminoAcids::amino_acids[amino_acid];
}

std::map<std::string, float> getAminoAcids() {
    return AminoAcids::amino_acids;
}
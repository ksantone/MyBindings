#include <map>
#include "AminoAcids.h"

std::map <std::string, float> AminoAcids::amino_acids;

/*std::map <float, std::string> amino_acids_backwards = { { 71.03711, "A"}, { 156.10111, "R" }, { 114.04293, "N" }, { 115.02694, "D" }, { 103.00919, "C" }, { 129.04259, "E" }, { 128.05858, "Q" },
        { 57.02146, "G" }, { 137.05891, "H" }, { 113.08406, "I/L" }, { 128.09496, "K" }, { 131.04049, "M" }, { 147.06841, "F" }, { 97.05276, "P" },
        { 87.03203, "S"}, { 101.04768, "T"}, { "W", 186.07931}, { "Y", 163.06333}, { "V", 99.06841} };*/

float getAminoAcidMass(std::string amino_acid) {
	return AminoAcids::amino_acids[amino_acid];
}

std::map<std::string, float> getAminoAcids() {
    return AminoAcids::amino_acids;
}
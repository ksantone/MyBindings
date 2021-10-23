#include <map>
#include <string>

#ifndef AminoAcids_HEADER
#define AminoAcids_HEADER

class AminoAcids {
public:
    static std::map <std::string, float> amino_acids;
    std::map<std::string, float> getAminoAcids();
    float getAminoAcidMass(std::string amino_acid);
    AminoAcids();
};

#endif
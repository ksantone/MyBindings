#include <string>

#ifndef Peak_HEADER
#define Peak_HEADER

class Peak {
public:
	float mass, charge, intensity, rt;
	Peak(float mass, int charge, float intensity, float rt);/* {
		this->mz = mz;
		this->intensity = intensity;
		this->rt = rt;
	};*/
};

#endif
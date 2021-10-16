#include <string>

#ifndef Peak_HEADER
#define Peak_HEADER

class Peak {
public:
	float mz, intensity, rt;
	Peak(float mz, float intensity, float rt);/* {
		this->mz = mz;
		this->intensity = intensity;
		this->rt = rt;
	};*/
};

#endif
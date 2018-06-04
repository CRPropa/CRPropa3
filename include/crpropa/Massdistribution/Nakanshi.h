#ifndef CRPROPA_NAKANSHI_H
#define CRPROPA_NAKANSHI_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <string>
#include <math.h>

namespace crpropa {

class Nakanshi {
/*
 @class Nakanshi
 @brief Modell for HI arXiv:astro-ph/0304338
	Modell for H2 arxiv:astro-ph/0610769
*/ 

private:
	bool isforHI;
	bool isforHII;
	bool isforH2;

public:

	Nakanshi();
	double getDensity(const Vector3d &position);
	double getHIDensity(const Vector3d &position);
	double getH2Density(const Vector3d &position);

	double getHIScaleheight(const Vector3d &position);
	double getHIPlanedensity(const Vector3d &position);

	double getH2Scaleheight(const Vector3d &position);
	double getH2Planedensity(const Vector3d &position);


	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	std::string getDiscription();

};

} //namespace crpropa

#endif //CRPROPA_NAKANSHI_H

		


#ifndef CRPROPA_NAKANISHI_H
#define CRPROPA_NAKANISHI_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"


#include <math.h>
#include <sstream>

#include "kiss/logger.h"

namespace crpropa {

class Nakanishi: public Density{
/*
 @class Nakanishi
 @brief Modell for HI arXiv:astro-ph/0304338
	Modell for H2 arxiv:astro-ph/0610769
	fit of the models given in arXiv:1607.07886
*/ 

private:
	bool isforHI = true;
	bool isforHII = false;
	bool isforH2 = true;

public:

	double getDensity(const Vector3d &position)const;
	double getNucleonDensity(const Vector3d &position)const;
	double getHIDensity(const Vector3d &position)const;
	double getH2Density(const Vector3d &position)const;

	double getHIScaleheight(const Vector3d &position)const;
	double getHIPlanedensity(const Vector3d &position)const;

	double getH2Scaleheight(const Vector3d &position)const;
	double getH2Planedensity(const Vector3d &position)const;


	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void setisforHI(bool HI);
	void setisforH2(bool H2);
	
	std::string getDescription();
};

} //namespace crpropa

#endif //CRPROPA_NAKANISHI_H

		


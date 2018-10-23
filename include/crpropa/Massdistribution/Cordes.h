#ifndef CRPROPA_CORDES_H
#define CRPROPA_CORDES_H

#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <math.h>
#include <sstream>

#include "kiss/logger.h"

namespace crpropa {
/*
	@class Cordes
	@brief Density for HII 
	Cordes et al., 1991, Nature 353,737
	*/

class Cordes: public Density {

private:
	// DO NOT CHANGE model typ!
	bool isforHI = false;	 
	bool isforHII = true;
	bool isforH2 = false;	

public:

	double getDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getNucleonDensity(const Vector3d &position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	std::string getDescription();
	
};

} //namespace

#endif //CRPROPA_CORDES_H

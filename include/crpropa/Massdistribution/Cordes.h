#ifndef CRPROPA_CORDES_H
#define CRPROPA_CORDES_H

#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <math.h>

namespace crpropa {
/*
	@class Cordes
	@brief Density for HII 
	Cordes et al., 1991, Nature 353,737
	*/

class Cordes: public Density {

private:
	bool isforHI = false;	//model typ DO NOT CHANGE!
	bool isforHII = true;
	bool isforH2 = false;	

public:

	Cordes();
	double getDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
};

} //namespace

#endif //CRPROPA_CORDES_H

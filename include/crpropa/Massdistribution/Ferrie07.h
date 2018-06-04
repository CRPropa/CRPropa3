#ifndef CRPROPA_FERRIE07_H
#define CRPROPA_FERRIE07_H

#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <math.h>

namespace crpropa {
/**
 @class FerrieI
 @brief distribution for Ferrie 2007: 
	Außen: 	ApJ, 497, 759
	Innen:	arxiv:	astro-ph/0702532
*/

class Ferrie: public Density {

private:

	bool isforHI = true;		// standard for all kind of distribution
	bool isforHII = true;
	bool isforH2 = true;
	double Rsun = 8500*pc;
	
	Vector3d CMZTrafo(const Vector3d &position) const; 
	Vector3d DISKTrafo(const Vector3d &position) const;

public:

	Ferrie();
	
    double getDensity( const Vector3d &position)const;
	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;

	void setisforHI(bool HI);
	void setisforHII(bool HII);
	void setisforH2(bool H2);
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
};

}//namespace crpropa

#endif //CRPROPA_FERRIE07_H



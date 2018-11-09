#ifndef CRPROPA_FERRIERE_H
#define CRPROPA_FERRIERE_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <sstream>

#include "kiss/logger.h"

namespace crpropa {
/**
 @class Ferriere 
 @brief distribution of hydrogen in the Milky Way 
 * Here in model Ferriere 2007
 * seperated in 2 regions (inner, outer). The border is for R=3 kpc in galactocentric radius. 
 * model is discribed in 
outer: ApJ, 497, 759
inner:	arxiv:	astro-ph/0702532
*/
class Ferriere: public Density {

private:
	// standard for all types of distribution
	bool isforHI = true;		
	bool isforHII = true;
	bool isforH2 = true;
	double Rsun = 8500*pc;	// distance sun-galactic center
	
public:
	Vector3d CMZTrafo(const Vector3d &position) const; // coordinate trafo for the CentralMolecularZone Region
	Vector3d DISKTrafo(const Vector3d &position) const; // coordinate trafo for the region of the disk in galactic center

	double getDensity(const Vector3d &position) const;
	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
	double getNucleonDensity(const Vector3d &position) const;

	void setIsForHI(bool HI);
	void setIsForHII(bool HII);
	void setIsForH2(bool H2);
	
	bool getIsForHI();
	bool getIsForHII();
	bool getIsForH2();
	
	std::string getDescription();
};

}//namespace crpropa

#endif //CRPROPA_FERRIERE_H



#ifndef CRPROPA_CONSTANTDENSITY_H
#define CRPROPA_CONSTANTDENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <math.h>
#include <sstream>

#include "kiss/logger.h"


namespace crpropa {

/*
@class ConstantDensity 
@brief Density module for Constant densitys in HI, HII and H2 component. 
*/
class ConstantDensity: public Density {

private:
	// default mode: all density types set to 0 and no activ component
	double HIdensitynumber  = 0;	
	double HIIdensitynumber = 0;
	double H2densitynumber  = 0;
	
	bool isHI = false;
	bool isHII = false;
	bool isH2 = false;

public:
	ConstantDensity(double HI, double HII, double H2);
	
	double getDensity(const Vector3d &position) const;
	
	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
	double getNucleonDensity(const Vector3d &position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void setHI(bool activate, double densitynumber);
	void setHI(bool activate);	//change type status and keep densitynumber as it is
	void setHI(double densitynumber);//change densitynumber and keep type status as it is
	
	void setHII(bool activate, double densitynumber);
	void setHII(bool activate);
	void setHII(double densitynumber);
	
	void setH2(bool activate, double densitynumber);
	void setH2(bool activate);
	void setH2(double densitynumber);
	
	std::string getDescription();
};

} //namespace

#endif //CRPROPA_CONSTANTDENSITY_H



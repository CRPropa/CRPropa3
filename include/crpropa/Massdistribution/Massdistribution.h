#ifndef CRPROPA_MASSDISTRIBUTION_H
#define CRPROPA_MASSDISTRIBUTION_H


#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <string>
#include <vector>


namespace crpropa {

class Massdistribution: public Density {
	std::vector<ref_ptr<Density> >distributionList;		// Vector with information on distribuition for HI (0); HII (1); H2 (2)
bool isforHI=false;
bool isforHII=false;
bool isforH2=false;

public:
	double getDensity(const Vector3d &position) const;
	void add(Density &density);	

};
	
} //namespace crpropa

#endif //CRPROPA_MASSDISTRIBUTION_H



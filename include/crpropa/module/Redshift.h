#ifndef CRPROPA_REDSHIFT_H
#define CRPROPA_REDSHIFT_H

#include "crpropa/Module.h"

namespace crpropa {

/**
 @class Redshift
 @brief Updates redshift and applies adiabatic energy loss according to the travelled distance.
 */
class Redshift: public Module {
public:
	Redshift();
	void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_REDSHIFT_H

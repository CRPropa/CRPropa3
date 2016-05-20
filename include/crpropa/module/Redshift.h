#ifndef CRPROPA_REDSHIFT_H
#define CRPROPA_REDSHIFT_H

#include "crpropa/Module.h"

namespace crpropa {

/**
 @class Redshift
 @brief Updates redshift and applies adiabatic energy loss according to the traveled distance.
 */
class Redshift: public Module {
public:
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class FutureRedshift
 @brief Updates redshift and applies adiabatic energy loss according to the traveled distance. Extends to negative redshift values to allow for symmetric time windows around z=0.
 */
class FutureRedshift: public Module {
public:
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_REDSHIFT_H

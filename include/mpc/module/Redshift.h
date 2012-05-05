#ifndef REDSHIFT_H_
#define REDSHIFT_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class SimpleRedshift
 @brief Simple redshift calculation and adiabatic energy loss

 Uses the relation c*z = H(0)*D to approximately calculate small redshifts as function of distance to a given center.
 */
class SimpleRedshift {
private:
	Vector3d center;

public:
	SimpleRedshift(Vector3d center);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* REDSHIFT_H_ */

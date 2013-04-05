#ifndef MPC_REDSHIFT_H_
#define MPC_REDSHIFT_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class Redshift
 @brief Updates redshift and applies adiabatic energy loss according to the travelled distance.
 */
class Redshift: public Module {
public:
	Redshift();
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif // MPC_REDSHIFT_H_

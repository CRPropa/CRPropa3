#ifndef MPC_DEFLECTION_H_
#define MPC_DEFLECTION_H_

#include "mpc/Module.h"
#include "mpc/magneticField/MagneticField.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/PhasePoint.h"

namespace mpc {

/**
 @class DeflectionCK
 @brief Propagation through magnetic fields using the Cash-Karp integrator.

 This module solves the equations of motion of a charged particle when propagating through a magnetic field.\n
 It uses the Runge-Kutta integration method with Cash-Karp coefficients.\n
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 */
class DeflectionCK: public Module {
public:
	ref_ptr<MagneticField> field;
	ExplicitRungeKutta<PhasePoint> erk;
	double tolerance;
	double minStep;

	DeflectionCK(ref_ptr<MagneticField> field, double tolerance = 1e-4,	double minimumStep = 0.1 * kpc);
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* MPC_DEFLECTION_H_ */


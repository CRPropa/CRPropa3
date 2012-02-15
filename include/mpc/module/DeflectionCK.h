#ifndef DEFLECTION_H_
#define DEFLECTION_H_

#include "mpc/Module.h"
#include "mpc/magneticField/magneticField.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/PhasePoint.h"

namespace mpc {

/**
 @class DeflectionCK
 @brief Propagation through magnetic fields using the Cash-Karp integrator.

 This module solves the equations of motion of a charged particle when propagating through a magnetic field.\n
 It uses the Runge-Kutta integration method with Cash-Karp coefficients.\n
 The step size control tries to keep the error close to, but smaller than the designated error.
 */
class DeflectionCK: public Module {
public:
	MagneticField *field;
	ExplicitRungeKutta<PhasePoint> erk;
	enum ControlType {
		NoStepSizeControl, WorstOffender, RMS
	};
	ControlType controlType;
	double tolerance;

	DeflectionCK(MagneticField *field, ControlType controlType,
			double tolerance);
	~DeflectionCK();
	std::string getDescription() const;
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
};

} // namespace mpc

#endif /* DEFLECTION_H_ */


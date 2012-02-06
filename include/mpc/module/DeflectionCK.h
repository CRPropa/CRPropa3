#ifndef DEFLECTION_H_
#define DEFLECTION_H_

#include "mpc/Module.h"
#include "mpc/magneticField/magneticField.hpp"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/PhasePoint.h"

namespace mpc {

/**
 * @class Deflection
 * @brief Magnetic deflection in 3D using a the Cash-Karp Runge-Kutta method
 * propagates the particle by a step particle.getNextStep() or smaller.
 * The step size control tries to keep the error close to, but smaller than the maxError
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

} /* namespace mpc */

#endif /* DEFLECTION_H_ */


#ifndef CRPROPA_DEFLECTION_H
#define CRPROPA_DEFLECTION_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/ExplicitRungeKutta.h"
#include "crpropa/PhasePoint.h"

namespace crpropa {

/**
 @class DeflectionCK
 @brief Propagation through magnetic fields using the Cash-Karp integrator.

 This module solves the equations of motion of a relativistic charged particle when propagating through a magnetic field.\n
 It uses the Runge-Kutta integration method with Cash-Karp coefficients.\n
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 Additionally a minimum and maximum size for the steps can be set.
 For neutral particles a rectalinear propagation is applied and a next step of the maximum step size proposed.
 */
class DeflectionCK: public Module {
private:
	ref_ptr<MagneticField> field;
	ExplicitRungeKutta<PhasePoint> erk;
	double tolerance; /*< target relative error of the numerical integration */
	double minStep; /*< minimum step size of the propagation */
	double maxStep; /*< maximum step size of the propagation */

public:
	DeflectionCK(ref_ptr<MagneticField> field = NULL, double tolerance = 1e-4,
			double minStep = 0.1 * kpc, double maxStep = 4000 * Mpc);
	void process(Candidate *candidate) const;
	void setField(ref_ptr<MagneticField> field);
	void setTolerance(double tolerance);
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);
	double getTolerance() const;
	double getMinimumStep() const;
	double getMaximumStep() const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_DEFLECTION_H


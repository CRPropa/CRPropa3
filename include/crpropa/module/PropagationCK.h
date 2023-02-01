#ifndef CRPROPA_PROPAGATIONCK_H
#define CRPROPA_PROPAGATIONCK_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"
#include "kiss/logger.h"

namespace crpropa {
/**
 * \addtogroup Propagation 
 * @{
 */

/**
 @class PropagationCK
 @brief Rectilinear propagation through magnetic fields using the Cash-Karp method.

 This module solves the equations of motion of a relativistic charged particle when propagating through a magnetic field.\n
 It uses the Runge-Kutta integration method with Cash-Karp coefficients.\n
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 Additionally a minimum and maximum size for the steps can be set.
 For neutral particles a rectilinear propagation is applied and a next step of the maximum step size proposed.
 */
class PropagationCK: public Module {
public:
	class Y {
	public:
		Vector3d x, u; /*< phase-point: position and direction */

		Y() {
		}

		Y(const Vector3d &x, const Vector3d &u) :
				x(x), u(u) {
		}

		Y(double f) :
				x(Vector3d(f, f, f)), u(Vector3d(f, f, f)) {
		}

		Y operator *(double f) const {
			return Y(x * f, u * f);
		}

		Y &operator +=(const Y &y) {
			x += y.x;
			u += y.u;
			return *this;
		}
	};

private:
	std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	ref_ptr<MagneticField> field;
	double tolerance; /*< target relative error of the numerical integration */
	double minStep; /*< minimum step size of the propagation */
	double maxStep; /*< maximum step size of the propagation */

public:
	/** Constructor for the adaptive Kash Carp.
	 * @param field
	 * @param tolerance	 tolerance is criterion for step adjustment. Step adjustment takes place only if minStep < maxStep
	 * @param minStep	   minStep/c_light is the minimum integration time step
	 * @param maxStep	   maxStep/c_light is the maximum integration time step. 
	 */
    PropagationCK(ref_ptr<MagneticField> field = NULL, double tolerance = 1e-4,
			double minStep = (0.1 * kpc), double maxStep = (1 * Gpc));

	void process(Candidate *candidate) const;

	// derivative of phase point, dY/dt = d/dt(x, u) = (v, du/dt)
	// du/dt = q*c^2/E * (u x B)
	Y dYdt(const Y &y, ParticleState &p, double z) const;

	void tryStep(const Y &y, Y &out, Y &error, double t,
			ParticleState &p, double z) const;

	void setField(ref_ptr<MagneticField> field);
	void setTolerance(double tolerance);
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);

	 /** get functions for the parameters of the class PropagationCK, similar to the set functions */
	ref_ptr<MagneticField> getField() const;
	
	/** get magnetic field vector at current candidate position
	 * @param pos   current position of the candidate
	 * @param z	 current redshift is needed to calculate the magnetic field
	 * @return	  magnetic field vector at the position pos */
	Vector3d getFieldAtPosition(Vector3d pos, double z) const;

	double getTolerance() const;
	double getMinimumStep() const;
	double getMaximumStep() const;
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PROPAGATIONCK_H

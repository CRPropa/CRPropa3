#ifndef CRPROPA_DIFFUSIONSDE_H
#define CRPROPA_DIFFUSIONSDE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/advectionField/AdvectionField.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include "kiss/logger.h"

namespace crpropa {
/**
 * \addtogroup Propagation
 * @{
 */

/**
 @class DiffusionSDE
 @brief Propagates candidates as pseudo(!)-particles.
 The time integration of SDEs is used to solve the transport equation.
 * Here an Euler-Mayurama integration scheme is used. The diffusion tensor
 * can be anisotropic with respect to the magnetic field line coordinates.
 * The integration of field lines is done via the CK-algorithm.
 */


class DiffusionSDE : public Module{

private:
	    ref_ptr<MagneticField> magneticField;
	    ref_ptr<AdvectionField> advectionField;
	    double minStep; // minStep/c_light is the minimum integration timestep
	    double maxStep; // maxStep/c_light is the maximum integration timestep
	    double tolerance; // tolerance is criterion for step adjustment. Step adjustment takes place when the tangential vector of the magnetic field line is calculated.
	    double epsilon; // ratio of parallel and perpendicular diffusion coefficient D_par = epsilon*D_perp
	    double alpha; // power law index of the energy dependent diffusion coefficient: D\propto E^alpha
	    double scale; // scaling factor for the diffusion coefficient D = scale*D_0

public:
	/** Constructor
	 @param magneticField	the magnetic field to be used 
	 @param tolerance		Tolerance is criterion for step adjustment. Step adjustment takes place when the  tangential vector of the magnetic field line is calculated.
	 @param minStep			minStep/c_light is the minimum integration time step
	 @param maxStep			maxStep/c_light is the maximum integration time step
	 @param epsilon			Ratio of parallel and perpendicular diffusion coefficient D_par = epsilon*D_perp
	 */
	DiffusionSDE(ref_ptr<crpropa::MagneticField> magneticField, double tolerance = 1e-4, double minStep = 10 * pc, double maxStep = 1 * kpc, double epsilon = 0.1);
	/** Constructor
	 @param magneticField	the magnetic field to be used 
	 @param advectionField	object containing advection field
	 @param tolerance		Tolerance is criterion for step adjustment. Step adjustment takes place when the  tangential vector of the magnetic field line is calculated.
	 @param minStep			minStep/c_light is the minimum integration time step
	 @param maxStep			maxStep/c_light is the maximum integration time step
	 @param epsilon			Ratio of parallel and perpendicular diffusion coefficient D_par = epsilon*D_perp
	 */
	DiffusionSDE(ref_ptr<crpropa::MagneticField> magneticField, ref_ptr<crpropa::AdvectionField> advectionField, double tolerance = 1e-4, double minStep = 10 * pc, double maxStep = 1 * kpc, double epsilon = 0.1);

	void process(crpropa::Candidate *candidate) const;

	void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &PosErr, double z, double propStep ) const;
	void driftStep(const Vector3d &Pos, Vector3d &LinProp, double h) const;
	void calculateBTensor(double rig, double BTen[], Vector3d pos, Vector3d dir, double z) const;

	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);
	void setTolerance(double tolerance);
	void setEpsilon(double kappa);
	void setAlpha(double alpha);
	void setScale(double Scale);
	void setMagneticField(ref_ptr<crpropa::MagneticField> magneticField);
	void setAdvectionField(ref_ptr<crpropa::AdvectionField> advectionField);

	double getMinimumStep() const;
	double getMaximumStep() const;
	double getTolerance() const;
	double getEpsilon() const;
	double getAlpha() const;
	double getScale() const;
	std::string getDescription() const;
  
  ref_ptr<MagneticField> getMagneticField() const;
	/** get magnetic field vector at current candidate position
	 @param pos   current position of the candidate
	 @param z	 current redshift is needed to calculate the magnetic field
	 @return	  magnetic field vector at the position pos */
	Vector3d getMagneticFieldAtPosition(Vector3d pos, double z) const;
	ref_ptr<AdvectionField> getAdvectionField() const;
	/** get advection field vector at current candidate position
	 @param pos   current position of the candidate
	 @return	  magnetic field vector at the position pos */
	Vector3d getAdvectionFieldAtPosition(Vector3d pos) const;

};
/** @}*/

} //namespace crpropa

#endif // CRPROPA_DIFFUSIONSDE_H

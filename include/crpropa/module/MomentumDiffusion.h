#ifndef CRPROPA_MOMENTUMDIFFUSION_H
#define CRPROPA_MOMENTUMDIFFUSION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include "kiss/logger.h"

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */



/**
 @class MomentumDiffusion
 @brief Changes the energy/absolute momentum of the pseudo(!)-particles
 * Here an Euler-Maruyama integration scheme is used. The momentum diffusion
 * scales with the parallel spatial diffusion coefficient as Dpp*Dxx \propto p^2 * v_A^2
 * Here, v_A is the Alven speed of the wave turbulence. 
 */

class MomentumDiffusion: public Module {

private:
	double alpha; // power law index of the energy dependent diffusion coefficient: D\propto E^alpha
	double scale; // scaling factor for the diffusion coefficient D = scale*D_0
	double vA; //Alvén speed
	double limit; // fraction of allowed mean energy change per step

public:
	/** Constructor
	@param alpha 	spectral index of the diffusion. Should match alpha in DiffusionSDE
	@param Scale 	Scaling of the diffusion coefficient. Should match alpha in DiffusionSDE
	@param aV	Alven speed of the wave turbulence
	*/

	MomentumDiffusion(double alpha, double scale);
	MomentumDiffusion();

	void process(Candidate *candidate) const;
	double calculateBScalar(double rig, double p) const;
	double calculateAScalar(double rig,  double p) const;

	void setAlpha(double alpha);
	void setScale(double Scale);
	void setAlvenSpeed(double aV);
	void setLimit(double limit);

	double getAlpha() const;
	double getScale() const;
	double getAlvenSpeed() const;
	double getLimit() const;

	std::string getDescription() const;

	
};

/** @}*/

}; //end namespace crpropa

#endif // CRPROPA_MOMENTUMDIFFUSION_H



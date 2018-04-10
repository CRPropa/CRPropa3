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


class MomentumDiffusion: public Module {

private:
	double alpha; // power law index of the energy dependent diffusion coefficient: D\propto E^alpha
	double scale; // scaling factor for the diffusion coefficient D = scale*D_0

public:
	MomentumDiffusion(double alpha, double scale);
	MomentumDiffusion();

	void process(Candidate *candidate) const;
	void calculateBScalar(double rig, double BScal) const;
	void calculateAScalar(double rig,  double p, int id, double AScal) const;

	void setAlpha(double alpha);
	void setScale(double Scale);

	double getAlpha() const;
	double getScale() const;

	
};



} //end namespace crpropa

#endif // CRPROPA_MOMENTUMDIFFUSION_H



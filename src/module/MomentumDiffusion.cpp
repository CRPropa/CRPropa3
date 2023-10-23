#include "crpropa/module/MomentumDiffusion.h"

using namespace crpropa;

ConstantMomentumDiffusion::ConstantMomentumDiffusion(double Dpp) {
	setLimit(0.1);
	setDpp(Dpp);
}

ConstantMomentumDiffusion::ConstantMomentumDiffusion(double Dpp, double limit) {
	setLimit(limit);
	setDpp(Dpp);
}

void ConstantMomentumDiffusion::process(Candidate *c) const {
	double rig = c->current.getRigidity();
	if (std::isinf(rig)) {
		return; // Only charged particles
	}
	
	double p = c->current.getEnergy() / c_light; // Note we use E=p/c (relativistic limit)
	double dt = c->getCurrentStep() / c_light;
	
	double eta =  Random::instance().randNorm();
	double domega = eta * sqrt(dt);
	
	double AScal = calculateAScalar(p);
	double BScal = calculateBScalar();

	double dp = AScal * dt + BScal * domega;
	c->current.setEnergy((p + dp) * c_light);
	
	c->limitNextStep(limit * p / AScal * c_light);
}

double ConstantMomentumDiffusion::calculateAScalar(double p) const {
	double a = + 2. / p * Dpp;
	return a; 
}

double ConstantMomentumDiffusion::calculateBScalar() const {
	double b = sqrt(2 * Dpp);
	return b;
}

void ConstantMomentumDiffusion::setDpp(double d) {
	if (d < 0 )
		throw std::runtime_error(
				"ConstantMomentumDiffusion: Dpp must be non-negative");
	Dpp = d;
}

void ConstantMomentumDiffusion::setLimit(double l) {
	limit = l;
}

double ConstantMomentumDiffusion::getDpp() const {
	return Dpp;
}

double ConstantMomentumDiffusion::getLimit() const {
	return limit;
}

std::string ConstantMomentumDiffusion::getDescription() const {
	std::stringstream s;
	s << "limit: " << limit << "\n";
	s << "Dpp: " << Dpp / (meter * meter / second) << " m^2/s";

	return s.str();
}

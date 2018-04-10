#include "crpropa/module/MomentumDiffusion.h"


using namespace crpropa;


MomentumDiffusion::MomentumDiffusion() : 
scale(1.)
{
  	setAlpha(1./2.);
	}

MomentumDiffusion::MomentumDiffusion(double alpha, double scale) :
scale(scale)
{
	setAlpha(alpha);
  	}

void MomentumDiffusion::process(Candidate *c) const {

	double p = c->current.getEnergy(); // Note we use E=p/c (relativistic limit)
	double rig = p / c->current.getCharge();
	int id = c->current.getId();
	double dt = c->getCurrentStep() / c_light;

	double eta =  Random::instance().randNorm();
	double domega = eta * pow(dt, 0.5);
	
	double AScal = 0.;
	double BScal = 0.;

	calculateAScalar(rig, p, id, AScal);
	calculateBScalar(rig, BScal);

	double dp = AScal * dt + BScal * domega;

	c->current.setEnergy(p + dp);

	//c->limitNextStep(limit * E / fabs(dEdt) *c_light);

}




void MomentumDiffusion::calculateBScalar(double rig, double BScal) const {

    double Dpp = scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha);
    BScal = pow( 2  * Dpp, 0.5);   
    return;

}

void MomentumDiffusion::calculateAScalar(double rig, double p, int id, double AScal) const {
	
	double Dpp = scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha);
    	double partialDpp = alpha * scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha-1);
	int sign = (id < 0) ? -1 : 1; // Different behavior for time forward (particles) and time backward (anti-particles) propagation  
    	AScal = partialDpp - sign * 2. / p * Dpp;   
   	return;

}


void MomentumDiffusion::setAlpha(double a) {
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"MomentumDiffusion: alpha not in range 0-2");
	alpha = a;
}

void MomentumDiffusion::setScale(double s) {
	if (s < 0)
		throw std::runtime_error(
				"MomentumDiffusion: Scale error: Scale < 0");
	scale = s;
}


double MomentumDiffusion::getAlpha() const {
	return alpha;
}

double MomentumDiffusion::getScale() const {
	return scale;
}



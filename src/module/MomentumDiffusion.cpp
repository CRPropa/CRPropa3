#include "crpropa/module/MomentumDiffusion.h"


using namespace crpropa;


MomentumDiffusion::MomentumDiffusion() : 
limit(0.1)
{
	setScale(1.);  	
	setAlpha(1./3.);
	setAlvenSpeed(40 * km / second);
	}

MomentumDiffusion::MomentumDiffusion(double alpha, double scale) :
limit(0.1)
{
	setScale(scale);
	setAlpha(alpha);
	setAlvenSpeed(40 * km / second);
  	}

void MomentumDiffusion::process(Candidate *c) const {

	double p = c->current.getEnergy()/c_light; // Note we use E=p/c (relativistic limit)
	double rig = p*c_light / c->current.getCharge();
	double dt = c->getCurrentStep() / c_light;
	
	std::cout <<dt<<"\n";
	double eta =  Random::instance().randNorm();
	double domega = eta * sqrt(dt);
	

	double AScal = calculateAScalar(rig, p);
	double BScal = calculateBScalar(rig, p);

	double dp = AScal * dt + BScal * domega;
	std::cout <<dp<<"\n";
	c->current.setEnergy((p + dp)*c_light);

	//c->limitNextStep(limit * p / ((AScal + BScal/sqrt(dt)) * c_light));
	c->limitNextStep(limit * p / ((AScal + BScal/sqrt(dt)) * c_light)); //Check for the factor c_light

}




double MomentumDiffusion::calculateBScalar(double rig, double p) const{

    	double Dxx = scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha);
	double Dpp = ( 4*vA*vA*p*p ) / ( 3*alpha*(4-alpha*alpha)*(4-alpha) ) / Dxx; // Astroparticle Physics: Theory and Phenomenology, G. Sigl, Atlantis Press (2017); Eq. (7.34) 
    	double BScal = sqrt( 2  * Dpp);   
    	return BScal;

}

// What is the physical interpretetation of this term? 7/27/19 LM
double MomentumDiffusion::calculateAScalar(double rig, double p) const {
	
    	double Dxx = scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha);
	double Dpp = ( 4*vA*vA*p*p ) / ( 3*alpha*(4-alpha*alpha)*(4-alpha) ) / Dxx;
    	double partialDpp = (2 - alpha) / p * Dpp; //check the sign: Should be correct 7/27/19 LM

    	double AScal = partialDpp -2. / p * Dpp; //=-alpha / p * Dpp
   	return AScal;

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

void MomentumDiffusion::setAlvenSpeed(double v) {
	if (v < 0)
		throw std::runtime_error(
				"MomentumDiffusion: Alvenspeed error: v_A < 0");
	if (v > c_light)
		throw std::runtime_error(
				"MomentumDiffusion: Alvenspeed error: v_A > c");
	vA = v;
}

void MomentumDiffusion::setLimit(double l) {
	if (l < 0)
		throw std::runtime_error(
				"MomentumDiffusion: Limit error: limit < 0");
	if (l > 1)
		throw std::runtime_error(
				"MomentumDiffusion: Limit error: limit > 1");
	limit = l;
}



double MomentumDiffusion::getAlpha() const {
	return alpha;
}

double MomentumDiffusion::getScale() const {
	return scale;
}

double MomentumDiffusion::getAlvenSpeed() const {
	return vA;
}

double MomentumDiffusion::getLimit() const {
	return limit;
}

std::string MomentumDiffusion::getDescription() const {
	std::stringstream s;
	s << "alpha: " << alpha << " , ";
	s << "Alven speed v_A: " << vA  << " km/s, ";
	s << "Scale: " << scale << "\n";

	return s.str();
}

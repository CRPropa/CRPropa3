// TO DO 12/12/19 LM
// 1.   Look into the limitNextStep function:
//      It is quite complicated calculation and it is not clear if it will help
//      Testing for this is required.
// 2.a) Write a test suite for the unit tests
// 	 b)	Write physics based tests
// 3. Provide simple examples


#include "crpropa/module/MomentumDiffusion.h"

using namespace crpropa;

ConstantMomentumDiffusion::ConstantMomentumDiffusion(double Dpp) {
	setDpp(Dpp);
}

void ConstantMomentumDiffusion::process(Candidate *c) const {
	double rig = c->current.getRigidity();
	if (std::isinf(rig)) {
		return; // Only charged particles
	}
	
	double p = c->current.getEnergy()/c_light; // Note we use E=p/c (relativistic limit)
	double dt = c->getCurrentStep() / c_light;
	
	std::cout <<dt<<"\n"; // Has to be deleted
	double eta =  Random::instance().randNorm();
	double domega = eta * sqrt(dt);
	
	double AScal = calculateAScalar(p);
	double BScal = calculateBScalar();

	double dp = AScal * dt + BScal * domega;
	std::cout <<dp<<"\n"; //Has to be deleted
	c->current.setEnergy((p + dp)*c_light);
	
	//Fast, but a little bit inconsitent, since the equation limit=\Delta p(step) / p
	//is not correctly solved.  
	//c->limitNextStep(limit * p / ((AScal + BScal/sqrt(dt)) * c_light)); //Check for the factor c_light
	
	//Solving the equation correctly
	//double c_tilde = limit*p;
	//double a = AScal;
	//double b = fabs(BScal)*0.8; //mean value of next diffusion step |eta|
	//c->limitNextStep((-b*sqrt(4*a*c_tilde+b*b)+4*a*c_tilde+b*b)/(2*a*a) * c_light);
}

double ConstantMomentumDiffusion::calculateAScalar(double p) const {
	double a = - 2./p * Dpp;
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

double ConstantMomentumDiffusion::getDpp() const {
	return Dpp;
}


PowerlawMomentumDiffusion::PowerlawMomentumDiffusion() : 
limit(0.1)
{
	setScale(1.);  	
	setAlpha(1./3.);
	setAlvenSpeed(40 * km / second);
	}

PowerlawMomentumDiffusion::PowerlawMomentumDiffusion(double alpha, double scale) :
limit(0.1)
{
	setScale(scale);
	setAlpha(alpha);
	setAlvenSpeed(40 * km / second);
  	}

void PowerlawMomentumDiffusion::process(Candidate *c) const {
	
	double rig = c->current.getRigidity();
	
	if (std::isinf(rig)) {
		return; // Only charged particles
	}
	
	double p = c->current.getEnergy()/c_light; // Note we use E=p/c (relativistic limit)
	double dt = c->getCurrentStep() / c_light;
	
	std::cout <<dt<<"\n"; // Has to be deleted
	double eta =  Random::instance().randNorm();
	double domega = eta * sqrt(dt);
	
	double Dpp = calculateDpp(rig, p);
	double AScal = calculateAScalar(rig, p, Dpp);
	double BScal = calculateBScalar(rig, p, Dpp);

	double dp = AScal * dt + BScal * domega;
	std::cout <<dp<<"\n"; //Has to be deleted
	c->current.setEnergy((p + dp)*c_light);
	
	//Fast, but a little bit inconsitent, since the equation limit=\Delta p(step) / p
	//is not correctly solved.  
	//c->limitNextStep(limit * p / ((AScal + BScal/sqrt(dt)) * c_light)); //Check for the factor c_light
	
	//Solving the equation correctly
	double c_tilde = limit*p;
	double a = AScal;
	double b = fabs(BScal)*0.8; //mean value of next diffusion step |eta|
	c->limitNextStep((-b*sqrt(4*a*c_tilde+b*b)+4*a*c_tilde+b*b)/(2*a*a) * c_light);
}

double PowerlawMomentumDiffusion::calculateDpp(double rig, double p) const {
	
	// The implementation of the spatial diffusion has to be the same as in DiffusionSDE
	double Dxx = scale * 6.1e24 * pow((std::abs(rig) / 4.0e9), alpha);
	
	// Astroparticle Physics: Theory and Phenomenology, G. Sigl, Atlantis Press (2017); Eq. (7.34)
	double Dpp = ( 4*vA*vA*p*p ) / ( 3*alpha*(4-alpha*alpha)*(4-alpha) ) / Dxx;
	 
	return Dpp;
}

double PowerlawMomentumDiffusion::calculateBScalar(double rig, double p, double Dpp) const{

    double BScal = sqrt( 2  * Dpp); 
      
    return BScal;
}

// What is the physical interpretetation of this term? 7/27/19 LM
double PowerlawMomentumDiffusion::calculateAScalar(double rig, double p, double Dpp) const {
	
    double partialDpp = (2 - alpha) / p * Dpp; //check the sign: Should be correct 7/27/19 LM
    double AScal = partialDpp -2. / p * Dpp; //=-alpha / p * Dpp
    
   	return AScal;
}

void PowerlawMomentumDiffusion::setAlpha(double a) {
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"MomentumDiffusion: alpha not in range 0-2");
	alpha = a;
}

void PowerlawMomentumDiffusion::setScale(double s) {
	if (s < 0)
		throw std::runtime_error(
				"MomentumDiffusion: Scale error: Scale < 0");
	scale = s;
}

void PowerlawMomentumDiffusion::setAlvenSpeed(double v) {
	if (v < 0)
		throw std::runtime_error(
				"MomentumDiffusion: Alvenspeed error: v_A < 0");
	if (v > c_light)
		throw std::runtime_error(
				"MomentumDiffusion: Alvenspeed error: v_A > c");
	vA = v;
}

void PowerlawMomentumDiffusion::setLimit(double l) {
	if (l < 0)
		throw std::runtime_error(
				"MomentumDiffusion: Limit error: limit < 0");
	if (l > 1)
		throw std::runtime_error(
				"MomentumDiffusion: Limit error: limit > 1");
	limit = l;
}

double PowerlawMomentumDiffusion::getAlpha() const {
	return alpha;
}

double PowerlawMomentumDiffusion::getScale() const {
	return scale;
}

double PowerlawMomentumDiffusion::getAlvenSpeed() const {
	return vA;
}

double PowerlawMomentumDiffusion::getLimit() const {
	return limit;
}

std::string PowerlawMomentumDiffusion::getDescription() const {
	std::stringstream s;
	s << "alpha: " << alpha << " , ";
	s << "Alven speed v_A: " << vA / (km/second)  << " km/s, ";
	s << "Scale: " << scale << "\n";

	return s.str();
}

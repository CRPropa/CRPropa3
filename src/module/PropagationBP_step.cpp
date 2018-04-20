#include "crpropa/module/PropagationBP_step.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {

PropagationBP_step::PropagationBP_step(ref_ptr<MagneticField> field,
		double minStep, double maxStep , double epsilon):minStep(0)
		 {
	setField(field);
	setMinStep(minStep);
	setMaxStep(maxStep);
	setEpsilon(epsilon);
}


void PropagationBP_step::process(Candidate *c) const {
	
	// update
	
	c->previous = c->current;	
	
	double step = clip(c->getNextStep() , minStep, maxStep);
	
	//get particle properties
	double q = c->current.getCharge();                    
	double m = c->current.getEnergy()/(c_light*c_light);  
	Vector3d x = c->current.getPosition();                
	Vector3d v = c->current.getDirection();
	double z = c->getRedshift();
	
	double EST;
	double rel_error;
	double alpha;
	
	int help = 1;
	
	do{
		
	//~ std::cout << "step/kpc =  " << step/kpc << std::endl;
	
	phasepoint p1 = testStep(q, m, x, v, z, step);  // 1 Schritt mit h
	
	phasepoint p_help = testStep(q, m, x, v, z, step/2);  // 2 Schritte mit h/2
	phasepoint p2 = testStep(q, m, p_help.x, p_help.v, z, step/2);
	
	//~ EST   = step * errorEstimation(p1.x , p2.x , step);
	EST   = errorEstimation(p1.x , p2.x , step);
	rel_error = EST/epsilon;
	
	alpha = 0.95 * pow(rel_error, - 1/3); // ( 3 = p + 1 ) // "-" wg rel_error = EST/TOL
	
	//~ std::cout << "error =  " << rel_error << ", alpha = " << alpha << ", NUM = " << help << std::endl;
	
	if(rel_error > 1){
		
		x = p1.x;
		v = p1.v;
		
		break;

		
	}else{
	
		step = clip(step * alpha, 0.1 * step, 5. * step );
		step = clip(step , minStep, maxStep);
		if(step == minStep){
			x = p1.x;
			v = p1.v;
			break;}
	}
	
	help ++;
	
	}while(1);
	
	
	c->current.setDirection(v);
	c->current.setPosition(x);
	
	c->setCurrentStep(step);     // Zeit -> Laenge
	c->setNextStep(step * alpha);
	
}

PropagationBP_step::phasepoint PropagationBP_step::testStep(double q, double m, Vector3d x, Vector3d v, double z, double step) const {
	
	PropagationBP_step::phasepoint p;
	double h = step / c_light;
	
	// half leap frog step in the position
	x += c_light * v * h /2. ;
	
	// get B field at particle position
	Vector3d B(0, 0, 0);
	try {
		B = field->getField(x, z);
	} catch (std::exception &e) {
		std::cerr << "PropagationBP: Exception in getField." << std::endl;
		std::cerr << e.what() << std::endl;
	}
	// Boris help vectors 
	Vector3d t = B * q/2/m * h;
	Vector3d s = t *2. /(1+t.dot(t));
	Vector3d v_help;
	
	// Boris Push
	v_help = v + v.cross(t);
	v = v + v_help.cross(s);
	
	// the other half leap frog step in the position
	x += c_light * v * h /2. ;
	
	p.x = x;
	p.v = v;
	
	return p;
	
}

double PropagationBP_step::errorEstimation(const Vector3d mu, const Vector3d muh, double h) const {

	Vector3d diff = (mu - muh);
	//~ std::cout << " diff = " << diff.getR() << std::endl;
	double S = diff.getR() / (h * (1 - 1/4 ) );    // 1/4 = (1/2)² = mu hoch p
	
	return S;
}


void PropagationBP_step::setField(ref_ptr<MagneticField> f) {
	field = f;
}

void PropagationBP_step::setMinStep(double min) {
	if (min < 0.)
		throw std::runtime_error("PropagationBP: Step < 0 ");
	minStep = min;
}

double PropagationBP_step::getMinStep() const {
	return minStep;
}

void PropagationBP_step::setMaxStep(double max) {
	if (max < 0.)
		throw std::runtime_error("PropagationBP: Step < 0 ");
	maxStep = max;
}

double PropagationBP_step::getMaxStep() const {
	return maxStep;
}

void PropagationBP_step::setEpsilon(double eps) {
	if (eps <= 0.)
		throw std::runtime_error("PropagationBP: maximum error delta < 0 ");
	epsilon = eps;
}

double PropagationBP_step::getEpsilon() const {
	return epsilon;
}

}

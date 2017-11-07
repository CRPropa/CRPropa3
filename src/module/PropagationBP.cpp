#include "crpropa/module/PropagationBP.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {

PropagationBP::PropagationBP(ref_ptr<MagneticField> field,
		double Step):Step(1* Mpc)
		 {
	setField(field);
	setStep(Step);
}


void PropagationBP::process(Candidate *c) const {
	
	// update
	c->previous = c->current;
	

	double step = getStep();
	c->setCurrentStep(step);
	step = step / c_light;
	
	//get particle properties
	double q = c->current.getCharge();                    
	double m = c->current.getEnergy()/(c_light*c_light);  
	Vector3d x = c->current.getPosition();                
	Vector3d v = c->current.getDirection();
	
	// half leap frog step in the position
	c->current.setPosition(x + c_light * v * step  ); 
	
	// get B field at particle position
	//Vector3d B(1 * nG, 0, 0);
	Vector3d B(0, 0, 0);
	//~ double z = 1; // <- Lukas fragen!
	double z = c->getRedshift();
	try {
		B = field->getField(x, z);
	} catch (std::exception &e) {
		std::cerr << "PropagationBP: Exception in getField." << std::endl;
		std::cerr << e.what() << std::endl;
	}
	// Boris help vectors 
	Vector3d t = B * q/2/m * step;
	Vector3d s = t *2. /(1+t.dot(t));
	Vector3d v_help;
	
	// Boris Push
	v_help = v + v.cross(t);
	v = v + v_help.cross(s);
	
	// full leap frog step in the velocity 
	c->current.setDirection(v);
	
	//~ // the other half leap frog step in the position
	c->current.setPosition(x +  c_light * v * step );
	
	c->setCurrentStep(step);
	c->setNextStep(step);
	
}


void PropagationBP::setField(ref_ptr<MagneticField> f) {
	field = f;
}

void PropagationBP::setStep(double min) {
	if (min < 0.)
		throw std::runtime_error("PropagationBP: Step < 0 ");
	Step = min;
}

double PropagationBP::getStep() const {
	return Step;
}

}

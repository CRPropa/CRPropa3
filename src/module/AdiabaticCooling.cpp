#include "crpropa/module/AdiabaticCooling.h"

namespace crpropa {

AdiabaticCooling::AdiabaticCooling(ref_ptr<AdvectionField> advectionField) :
	advectionField(advectionField) {
	setLimit(0.1);
}

AdiabaticCooling::AdiabaticCooling(ref_ptr<AdvectionField> advectionField, double limit) :
	advectionField(advectionField) {
	setLimit(limit);
}

void AdiabaticCooling::process(Candidate *c) const {

	Vector3d pos = c->current.getPosition();
	double E = c->current.getEnergy(); // Note we use E=p/c (relativistic limit)
	
	double Div = 0.;
	
	try {
		Div = advectionField->getDivergence(pos);
	} 
	catch (std::exception &e) {
		std::cerr << "AdiabaticCooling: Exception in getDivergence." << std::endl;
		std::cerr << e.what() << std::endl;
	}

	double dEdt = -E / 3. * Div; 	// cooling due to advection 
					// (see e.g. Kopp et al. Computer Physics Communication 183
					// (201) 530-542)
	double dt = c->getCurrentStep() / c_light;
	double dE = dEdt * dt;
	
	c->current.setEnergy(E + dE);
	c->limitNextStep(limit * E / dEdt *c_light);
}

void AdiabaticCooling::setLimit(double l) {
	limit = l;
	return;
}

double AdiabaticCooling::getLimit() const {
	return limit;
}
	
	



} // end namespace crpropa

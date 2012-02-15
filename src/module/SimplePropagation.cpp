#include "mpc/module/SimplePropagation.h"

namespace mpc {

SimplePropagation::SimplePropagation() :
		acceleration(5) {
}

SimplePropagation::SimplePropagation(double acceleration) :
		acceleration(acceleration) {
}

void SimplePropagation::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	double nextStep = candidate->getNextStep();
	Vector3 pos = candidate->current.getPosition();
	Vector3 dir = candidate->current.getDirection();
	candidate->current.setPosition(pos + dir * nextStep);
	candidate->setCurrentStep(nextStep);
	candidate->setNextStep(nextStep * acceleration);
}

std::string SimplePropagation::getDescription() const {
	return "Simple Propagation";
}

} // namespace mpc

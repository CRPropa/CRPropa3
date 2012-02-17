#include "mpc/module/SimplePropagation.h"

namespace mpc {

SimplePropagation::SimplePropagation() :
		acceleration(10), initialStep(10 * kpc) {
}

SimplePropagation::SimplePropagation(double acceleration, double initialStep) :
		acceleration(acceleration), initialStep(initialStep) {
}

void SimplePropagation::process(Candidate *candidate) {
	double nextStep = candidate->getNextStep();
	if (nextStep == 0)
		nextStep = initialStep;
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

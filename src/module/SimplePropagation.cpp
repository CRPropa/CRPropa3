#include "mpc/module/SimplePropagation.h"

namespace mpc {

SimplePropagation::SimplePropagation(double acceleration, double minimumStep) :
		acceleration(acceleration), minimumStep(minimumStep) {
}

void SimplePropagation::process(Candidate *candidate) const {
	double nextStep = std::max(minimumStep, candidate->getNextStep());
	Vector3d pos = candidate->current.getPosition();
	Vector3d dir = candidate->current.getDirection();
	candidate->current.setPosition(pos + dir * nextStep);
	candidate->setCurrentStep(nextStep);
	candidate->setNextStep(nextStep * acceleration);
}

std::string SimplePropagation::getDescription() const {
	return "Simple Propagation";
}

} // namespace mpc

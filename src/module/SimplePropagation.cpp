#include "mpc/module/SimplePropagation.h"

namespace mpc {

SimplePropagation::SimplePropagation(double minStep, double maxStep) :
		minStep(minStep), maxStep(maxStep) {
}

void SimplePropagation::setMinimumStep(double s) {
	minStep = s;
}

void SimplePropagation::setMaximumStep(double s) {
	maxStep = s;
}

void SimplePropagation::process(Candidate *candidate) const {
	// save the new previous particle state
	candidate->previous = candidate->current;

	double step = candidate->getNextStep();
	step = std::max(step, minStep);

	Vector3d pos = candidate->current.getPosition();
	Vector3d dir = candidate->current.getDirection();
	candidate->current.setPosition(pos + dir * step);
	candidate->setCurrentStep(step);
	candidate->setNextStep(maxStep);
}

std::string SimplePropagation::getDescription() const {
	std::stringstream s;
	s << "Simple Propagation: Step size = " << minStep / Mpc
			<< " - " << maxStep / Mpc << " Mpc";
	return s.str();
}

} // namespace mpc

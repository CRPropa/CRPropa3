#include "mpc/module/SimplePropagation.h"

namespace mpc {

SimplePropagation::SimplePropagation(double minStep, double maxStep) :
		minStep(minStep), maxStep(maxStep) {
}

void SimplePropagation::process(Candidate *candidate) const {
	double step = candidate->getNextStep();
	step = std::max(step, minStep);
	Vector3d pos = candidate->current.getPosition();
	Vector3d dir = candidate->current.getDirection();
	candidate->current.setPosition(pos + dir * step);
	candidate->setCurrentStep(step);
	candidate->setNextStep(maxStep);
}

std::string SimplePropagation::getDescription() const {
	return "Simple Propagation";
}

} // namespace mpc

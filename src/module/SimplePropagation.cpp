#include "mpc/module/SimplePropagation.h"

#include <sstream>
#include <stdexcept>

namespace mpc {

SimplePropagation::SimplePropagation(double minStep, double maxStep) :
		minStep(minStep), maxStep(maxStep) {
	if (minStep > maxStep)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
}

void SimplePropagation::process(Candidate *c) const {
	c->previous = c->current;

	double step = std::max(minStep, c->getNextStep());
	c->setCurrentStep(step);

	Vector3d pos = c->current.getPosition();
	Vector3d dir = c->current.getDirection();
	c->current.setPosition(pos + dir * step);

	c->setNextStep(maxStep);
}

void SimplePropagation::setMinimumStep(double step) {
	if (step > maxStep)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
	minStep = step;
}

void SimplePropagation::setMaximumStep(double step) {
	if (minStep > step)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
	maxStep = step;
}

double SimplePropagation::getMinimumStep() const {
	return minStep;
}

double SimplePropagation::getMaximumStep() const {
	return maxStep;
}

std::string SimplePropagation::getDescription() const {
	std::stringstream s;
	s << "Simple Propagation: Step size = " << minStep / Mpc
			<< " - " << maxStep / Mpc << " Mpc";
	return s.str();
}

} // namespace mpc

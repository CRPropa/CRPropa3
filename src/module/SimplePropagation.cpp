#include "crpropa/module/SimplePropagation.h"

#include <sstream>
#include <stdexcept>

namespace crpropa {

SimplePropagation::SimplePropagation(double minStep, double maxStep) :
		minStep(minStep/c_light), maxStep(maxStep/c_light) {
	if (minStep > maxStep)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
}

SimplePropagation::SimplePropagation(bool doTimeSteps, double minStep, double maxStep){
	if (doTimeSteps){
		this->minStep = minStep;
		this->maxStep = maxStep;
	} else {
		this->minStep = minStep/c_light;
		this->maxStep = maxStep/c_light;
	}

	if (minStep > maxStep)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
}

void SimplePropagation::process(Candidate *c) const {
	c->previous = c->current;

	double dt = clip(c->getNextStep(), minStep, maxStep);
	c->setCurrentStep(dt);
	Vector3d pos = c->current.getPosition();
	Vector3d dir = c->current.getDirection();
	double vel = c->getVelocity();
	c->current.setPosition(pos + dir * vel * dt);
	c->setNextStep(maxStep);
}

void SimplePropagation::setMinimumStep(double step) {
	if (step/c_light > maxStep)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
	minStep = step/c_light;
}

void SimplePropagation::setMaximumStep(double step) {
	if (minStep > step/c_light)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
	maxStep = step/c_light;
}

void SimplePropagation::setMinimumTimeStep(double step){
	if (step > maxStep)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
	minStep = step;
}

void SimplePropagation::setMaximumTimeStep(double step){
	if (minStep > step)
		throw std::runtime_error("SimplePropagation: minStep > maxStep");
	maxStep = step;
}

std::string SimplePropagation::getDescription() const {
	std::stringstream s;
	s << "SimplePropagation: Step size = " << minStep / kiloyear
			<< " - " << maxStep / kiloyear << " kiloyear";
	return s.str();
}

} // namespace crpropa

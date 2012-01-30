#include "mpc/Candidate.h"

namespace mpc {

Candidate::Candidate() :
		redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), status(Active) {
}

double Candidate::getRedshift() const {
	return redshift;
}

double Candidate::getTrajectoryLength() const {
	return trajectoryLength;
}

double Candidate::getCurrentStep() const {
	return currentStep;
}

double Candidate::getNextStep() const {
	return nextStep;
}

Candidate::Status Candidate::getStatus() const {
	return status;
}

void Candidate::setRedshift(double z) {
	redshift = z;
}

void Candidate::setTrajectoryLength(double a) {
	trajectoryLength = a;
}

void Candidate::setCurrentStep(double lstep) {
	currentStep = lstep;
	trajectoryLength += lstep;
}

void Candidate::setNextStep(double step) {
	nextStep = step;
}

void Candidate::limitNextStep(double step) {
	nextStep = std::min(nextStep, step);
}

void Candidate::setStatus(Status stat) {
	status = stat;
}

void Candidate::addSecondary(int id, double energy) {
	ParticleState p = current; // makes a copy, right?
	p.setId(id);
	p.setEnergy(energy);
	Candidate c;
	c.status = Active;
	c.initial = p;
	c.current = p;
	c.redshift = redshift;
	c.trajectoryLength = trajectoryLength;
	c.setNextStep(getCurrentStep());
	secondaries.push_back(c);
}

} // namespace mpc

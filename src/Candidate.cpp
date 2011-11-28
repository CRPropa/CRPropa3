#include "mpc/Candidate.h"

namespace mpc {

Candidate::Candidate() :
		age(0), currentStep(0), nextStep(0), status(Active) {
}

double Candidate::getTrajectoryLength() const {
	return age;
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

void Candidate::setTrajectoryLength(double a) {
	age = a;
}

void Candidate::setCurrentStep(double lstep) {
	currentStep = lstep;
	age += lstep;
}

void Candidate::setNextStep(double nstep) {
	nextStep = nstep;
}

void Candidate::setStatus(Status stat) {
	status = stat;
}

} // namespace mpc

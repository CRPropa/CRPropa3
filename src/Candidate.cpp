#include "mpc/Candidate.h"

namespace mpc {

Candidate::Candidate() :
		age(0), lastStepSize(0), nextStepSize(0), status(Active) {
}

double Candidate::getTrajectoryLength() const {
	return age;
}

double Candidate::getLastStep() const {
	return lastStepSize;
}

double Candidate::getNextStep() const {
	return nextStepSize;
}

Candidate::Status Candidate::getStatus() const {
	return status;
}

void Candidate::setTrajectoryLength(double a) {
	age = a;
}

void Candidate::setLastStepSize(double lstep) {
	lastStepSize = lstep;
	age += lstep;
}

void Candidate::setNextStepSize(double nstep) {
	nextStepSize = nstep;
}

void Candidate::setStatus(Status stat) {
	status = stat;
}

} // namespace mpc

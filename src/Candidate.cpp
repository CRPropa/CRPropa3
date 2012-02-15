#include "mpc/Candidate.h"

namespace mpc {

Candidate::Candidate() :
		redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), status(
				Active) {
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

bool Candidate::getInteractionState(std::string moduleName,
		InteractionState &state) {
	std::map<std::string, InteractionState>::iterator i =
			interactionStates.find(moduleName);
	if (i == interactionStates.end())
		return false;
	state = i->second;
	return true;
}

void Candidate::setInteractionState(std::string moduleName,
		InteractionState state) {
	interactionStates[moduleName] = state;
}

void Candidate::clearInteractionStates() {
	interactionStates.clear();
}

void addSecondary(std::vector<Candidate *> &secondaries, Candidate *parent,
		int id, double energy) {
	Candidate *c = new Candidate;
	c->setStatus(Candidate::Active);
	c->setRedshift(parent->getRedshift());
	c->setTrajectoryLength(parent->getTrajectoryLength());
	c->setNextStep(parent->getCurrentStep());

	ParticleState p = parent->current;
	p.setId(id);
	p.setEnergy(energy);
	c->initial = p;
	c->current = p;
	secondaries.push_back(c);
}

} // namespace mpc

#include "mpc/Candidate.h"

namespace mpc {

Candidate::Candidate() :
		redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), status(
				Active) {
}

Candidate::~Candidate() {
	clearSecondaries();
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

const std::map<std::string, InteractionState> Candidate::getInteractionStates() const {
	return interactionStates;
}

void Candidate::addSecondary(int id, double energy) {
	Candidate *secondary = new Candidate;
	secondary->setStatus(Candidate::Active);
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength);
	secondary->setNextStep(0);
	secondary->initial = initial;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
	secondaries.push_back(secondary);
}

void Candidate::clearSecondaries() {
	for (size_t i = 0; i < secondaries.size(); i++)
		delete secondaries[i];
	secondaries.clear();
}

} // namespace mpc

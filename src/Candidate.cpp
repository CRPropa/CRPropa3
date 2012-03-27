#include "mpc/Candidate.h"

namespace mpc {

Candidate::Candidate() :
		redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), active(
				true) {
}

Candidate::Candidate(const ParticleState &state) :
		current(state), initial(state), redshift(0), trajectoryLength(0), currentStep(
				0), nextStep(0), active(true) {
}

bool Candidate::isActive() const {
	return active;
}

void Candidate::setActive(const bool b) {
	active = b;
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

bool Candidate::getInteractionState(const std::string &moduleName,
		InteractionState &state) const {
	std::map<std::string, InteractionState>::const_iterator i =
			interactionStates.find(moduleName);
	if (i == interactionStates.end())
		return false;
	state = i->second;
	return true;
}

void Candidate::setInteractionState(const std::string &moduleName,
		const InteractionState &state) {
	interactionStates[moduleName] = state;
}

const std::map<std::string, InteractionState> Candidate::getInteractionStates() const {
	return interactionStates;
}

void Candidate::clearInteractionStates() {
	interactionStates.clear();
}

bool Candidate::getProperty(const std::string &name, std::string &value) const {
	std::map<std::string, std::string>::const_iterator i = properties.find(
			name);
	if (i == properties.end())
		return false;
	value = i->second;
	return true;
}

void Candidate::setProperty(const std::string &name, const std::string &value) {
	properties[name] = value;
}

bool Candidate::hasProperty(const std::string &name) const {
	std::map<std::string, std::string>::const_iterator i = properties.find(
			name);
	if (i == properties.end())
		return false;
	return true;
}
;

void Candidate::addSecondary(int id, double energy) {
	ref_ptr<Candidate> secondary = new Candidate;
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength);
	secondary->initial = initial;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
	secondaries.push_back(secondary);
}

void Candidate::clearSecondaries() {
	secondaries.clear();
}

} // namespace mpc

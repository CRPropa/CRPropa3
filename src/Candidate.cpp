#include "crpropa/Candidate.h"

namespace crpropa {

Candidate::Candidate() :
		redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), active(
				true) {
}

Candidate::Candidate(const ParticleState &state) :
		source(state), created(state), current(state), previous(state), redshift(
				0), trajectoryLength(0), currentStep(0), nextStep(0), active(
				true) {
}

bool Candidate::isActive() const {
	return active;
}

void Candidate::setActive(bool b) {
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

bool Candidate::getInteractionState(const std::string &name,
		InteractionState &state) const {
	InteractionStatesMap::const_iterator i = interactionStates.find(name);
	if (i == interactionStates.end())
		return false;
	state = i->second;
	return true;
}

const Candidate::InteractionStatesMap Candidate::getInteractionStates() const {
	return interactionStates;
}

void Candidate::setInteractionState(const std::string &name,
		const InteractionState &state) {
	interactionStates[name] = state;
}

void Candidate::removeInteractionState(const std::string &name) {
	InteractionStatesMap::iterator i = interactionStates.find(name);
	if (i == interactionStates.end())
		return;
	interactionStates.erase(i);
}

void Candidate::clearInteractionStates() {
	interactionStates.clear();
}

void Candidate::setProperty(const std::string &name, const std::string &value) {
	properties[name] = value;
}

bool Candidate::removeProperty(const std::string& name) {
	PropertyMap::iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	properties.erase(i);
	return true;
}

bool Candidate::getProperty(const std::string &name, std::string &value) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	value = i->second;
	return true;
}

bool Candidate::hasProperty(const std::string &name) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	return true;
}

void Candidate::addSecondary(int id, double energy) {
	ref_ptr<Candidate> secondary = new Candidate;
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength);
	secondary->source = source;
	secondary->created = current;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
	secondaries.push_back(secondary);
}

void Candidate::clearSecondaries() {
	secondaries.clear();
}

} // namespace crpropa

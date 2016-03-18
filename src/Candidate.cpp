#include "crpropa/Candidate.h"
#include "crpropa/Random.h"

namespace crpropa {

Candidate::Candidate(int id, double E, Vector3d pos, Vector3d dir, double z) :
		trajectoryLength(0), currentStep(0), nextStep(0), active(true) {
	ParticleState state(id, E, pos, dir);
	source = state;
	created = state;
	previous = state;
	current = state;
	setRedshift(z);
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

void Candidate::setProperty(const std::string &name, const std::string &value) {
	properties[name] = value;
}

bool Candidate::getProperty(const std::string &name, std::string &value) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	value = i->second;
	return true;
}

bool Candidate::removeProperty(const std::string& name) {
	PropertyMap::iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	properties.erase(i);
	return true;
}

bool Candidate::hasProperty(const std::string &name) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	return true;
}

void Candidate::addSecondary(Candidate *c) {
	secondaries.push_back(c);
}

void Candidate::addSecondary(int id, double energy) {
	ref_ptr<Candidate> secondary = new Candidate;
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength);
	secondary->source = source;
	secondary->previous = previous;
	secondary->created = current;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
	secondaries.push_back(secondary);
}

void Candidate::addSecondary(int id, double energy, Vector3d position) {
	ref_ptr<Candidate> secondary = new Candidate;
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength);
	secondary->source = source;
	secondary->previous = previous;
	secondary->created = current;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
  secondary->current.setPosition(position);
  secondary->created.setPosition(position);
	secondaries.push_back(secondary);
}

void Candidate::clearSecondaries() {
	secondaries.clear();
}

std::string Candidate::getDescription() const {
	std::stringstream ss;
	ss << "CosmicRay at z = " << getRedshift() << "\n";
	ss << "  source:  " << source.getDescription() << "\n";
	ss << "  current: " << current.getDescription();
	return ss.str();
}

ref_ptr<Candidate> Candidate::clone(bool recursive) const {
	ref_ptr<Candidate> cloned = new Candidate;
	cloned->source = source;
	cloned->created = created;
	cloned->current = current;
	cloned->previous = previous;

	cloned->properties = properties;
	cloned->active = active;
	cloned->redshift = redshift;
	cloned->trajectoryLength = trajectoryLength;
	cloned->currentStep = currentStep;
	cloned->nextStep = nextStep;
	if (recursive) {
		cloned->secondaries.reserve(secondaries.size());
		for (size_t i = 0; i < secondaries.size(); i++) {
			cloned->secondaries.push_back(secondaries[i]->clone(recursive));
		}
	}
	return cloned;
}

const Vector3d randomPositionInPropagationStep(Candidate *c) {
  Vector3d initial = c->current.getPosition();
  Vector3d final = c->previous.getPosition();
  Random &random = Random::instance();
  return initial + random.rand() * (final - initial);
}

} // namespace crpropa

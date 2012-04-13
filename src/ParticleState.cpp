#include "mpc/ParticleState.h"

namespace mpc {

ParticleState::ParticleState() :
		id(0), pmass(0), energy(0), position(0, 0, 0), direction(1, 0, 0) {

}

void ParticleState::setPosition(const Vector3 &pos) {
	position = pos;
}

const Vector3 &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3 &dir) {
	direction = dir / dir.mag();
}

const Vector3 &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = newEnergy;
}

double ParticleState::getEnergy() const {
	return energy;
}

void ParticleState::setId(int newId) {
	id = newId;
	pmass = getNucleusMass(id);
}

int ParticleState::getId() const {
	return id;
}

int ParticleState::getChargeNumber() const {
	return getChargeNumberFromNucleusId(id);
}

double ParticleState::getCharge() const {
	return getChargeNumberFromNucleusId(id) * eplus;
}

int ParticleState::getMassNumber() const {
	return getMassNumberFromNucleusId(id);
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getLorentzFactor() const {
	return energy / (this->getMass() * c_squared);
}

void ParticleState::setLorentzFactor(double gamma) {
	energy = gamma * this->getMass() * c_squared;
}

Vector3 ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3 ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

} // namespace mpc

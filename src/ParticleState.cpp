#include "mpc/ParticleState.h"

namespace mpc {

int getNucleusId(int a, int z) {
	return 1000000000 + z * 10000 + a * 10;
}

int getChargeNumberFromNucleusId(int id) {
	return (id - 1000000000) / 10000;
}

int getMassNumberFromNucleusId(int id) {
	return ((id - 1000000000) % 10000) / 10;
}

double ParticleState::getMass() const {
	return getChargeNumberFromNucleusId(id) * amu;
}

double ParticleState::getChargeNumber() const {
	return getChargeNumberFromNucleusId(id);
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getLorentzFactor() const {
	return energy / (this->getMass() * c_squared);
}

const Vector3 &ParticleState::getPosition() const {
	return position;
}

const Vector3 &ParticleState::getDirection() const {
	return direction;
}

Vector3 ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3 ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

int ParticleState::getId() const {
	return id;
}

void ParticleState::setDirection(const Vector3 &dir) {
	direction = dir / dir.mag();
}

void ParticleState::setPosition(const Vector3 &pos) {
	position = pos;
}

void ParticleState::setId(int newId) {
	id = newId;
}

void ParticleState::setMass(double newMass) {
	mass = newMass;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = newEnergy;
}

}

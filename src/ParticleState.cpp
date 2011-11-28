#include "mpc/ParticleState.h"

namespace mpc {

double ParticleState::getMass() const {
	return massNumber * amu;
}
double ParticleState::getChargeNumber() const {
	return chargeNumber;
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

void ParticleState::setDirection(const Vector3 &dir) {
	direction = dir / dir.mag();
}

void ParticleState::setPosition(const Vector3 &pos) {
	position = pos;
}

void ParticleState::setChargeNumber(size_t charge) {
	chargeNumber = charge;
}

void ParticleState::setMass(double newMass) {
	mass = newMass;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = newEnergy;
}

ParticleState::Type ParticleState::getType() {
	return type;
}

void ParticleState::setType(Type t) {
	type = t;
}

}

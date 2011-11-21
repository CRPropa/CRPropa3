#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"

namespace mpc {

class ParticleState {
public:
	enum Type {
		Gamma, Lepton, Hadron
	};

	union {
		size_t baryonNumber;
		size_t leptonNumber;
		size_t massNumber;
	};

	double getEnergy() const;
	void setEnergy(double newEnergy);

	const Vector3 &getPosition() const;
	void setPosition(const Vector3 &pos);

	const Vector3 &getDirection() const;
	void setDirection(const Vector3 &dir);

	double getChargeNumber() const;
	void setChargeNumber(size_t charge);

	double getMass() const;
	void setMass(double newMass);

	Type getType();
	void setType(Type t);

	// convenience
	double getLorentzFactor() const;
	Vector3 getVelocity() const;
	Vector3 getMomentum() const;

private:
	double energy;
	Vector3 position;
	Vector3 direction;
	size_t chargeNumber;
	double mass;
	Type type;

};

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
	direction = dir;
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

} // namespace mpc

#endif /* PARTICLE_H_ */

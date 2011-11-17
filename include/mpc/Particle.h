#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "mpc/ThreeVector.h"
#include "mpc/Units.h"

namespace mpc {

class Particle {
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

	const Hep3Vector &getPosition() const;
	void setPosition(const Hep3Vector &pos);

	const Hep3Vector &getDirection() const;
	void setDirection(const Hep3Vector &dir);

	double getChargeNumber() const;
	void setChargeNumber(size_t charge);

	double getMass() const;
	void setMass(double newMass);

	Type getType();
	void setType(Type t);

	// convenience
	double getLorentzFactor() const;
	Hep3Vector getVelocity() const;
	Hep3Vector getMomentum() const;

private:
	double energy;
	Hep3Vector position;
	Hep3Vector direction;
	size_t chargeNumber;
	double mass;
	Type type;

};

double Particle::getMass() const {
	return massNumber * amu;
}
double Particle::getChargeNumber() const {
	return chargeNumber;
}
double Particle::getEnergy() const {
	return energy;
}
double Particle::getLorentzFactor() const {
	return energy / (this->getMass() * c_squared);
}
const Hep3Vector &Particle::getPosition() const {
	return position;
}
const Hep3Vector &Particle::getDirection() const {
	return direction;
}
Hep3Vector Particle::getVelocity() const {
	return direction * c_light;
}
Hep3Vector Particle::getMomentum() const {
	return direction * (energy / c_light);
}

void Particle::setDirection(const Hep3Vector &dir) {
	direction = dir;
}

void Particle::setPosition(const Hep3Vector &pos) {
	position = pos;
}

void Particle::setChargeNumber(size_t charge) {
	chargeNumber = charge;
}

void Particle::setMass(double newMass) {
	mass = newMass;
}

void Particle::setEnergy(double newEnergy) {
	energy = newEnergy;
}

Particle::Type Particle::getType() {
	return type;
}

void Particle::setType(Type t) {
	type = t;
}

} // namespace mpc

#endif /* PARTICLE_H_ */

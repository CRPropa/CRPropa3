#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "ThreeVector.h"
#include "Units.h"
namespace mpc {
class Particle {
public:
	double energy;
	enum Type {
		Gamma, Lepton, Hadron
	};

	Type type;
	double mass;
	size_t chargeNumber;
	Hep3Vector position;
	Hep3Vector direction;
	union {
		size_t baryonNumber;
		size_t leptonNumber;
		size_t massNumber;
	};

	double getCharge() const {
		return chargeNumber * eplus;
	}
	double getMass() const {
		return massNumber * amu;
	}
	double getChargeNumber() const {
		return chargeNumber;
	}
	double getMassNumber() const {
		return massNumber;
	}
	double getEnergy() const {
		return energy;
	}
	double getEnergyEeV() const {
		return energy / EeV;
	}
	double getLorentzFactor() const {
		return energy / (this->getMass() * c_squared);
	}
	const Hep3Vector &getPosition() const {
		return position;
	}
	Hep3Vector getPositionMpc() const {
		return position / Mpc;
	}
	const Hep3Vector &getDirection() const {
		return direction;
	}
	Hep3Vector getVelocity() const {
		return direction * c_light;
	}
	Hep3Vector getMomentum() const {
		return direction * (energy / c_light);
	}

	void setDirection(const Hep3Vector &dir) {
		direction = dir;
	}

	void setPosition(const Hep3Vector &pos) {
		position = pos;
	}

	void setPositionMpc(const Hep3Vector &pos) {
		position = pos * Mpc;
	}

	void setChargeNumber(size_t charge) {
		chargeNumber = charge;
	}

	void setMass(double newMass) {
		mass = newMass;
	}

	void setEnergyEeV(double newEnergy) {

		energy = newEnergy * EeV;
	}
	void setEnergy(double newEnergy) {
		energy = newEnergy;
	}

	Type getType() {
		return type;
	}

	void setType(Type t) {
		type = t;
	}
};

class Candidate {
public:
	enum Status {
		Active = 0, Detected, AboveMaxTime, BelowEnergyThreshold, Decayed
	};

	Particle current;
	Particle initial;
	Candidate *parent;

	double age;
	double lastStep, nextStep;
	Status status;

	Candidate() :
			parent(0), age(0), lastStep(0), nextStep(0), status(Active) {

	}

	double getTrajectoryLength() const {
		return age;
	}
	double getTrajectoryLengthMpc() const {
		return age / Mpc;
	}
	double getLastStep() const {
		return lastStep;
	}
	double getLastStepMpc() const {
		return lastStep / Mpc;
	}
	double getNextStep() const {
		return nextStep;
	}
	double getNextStepMpc() const {
		return nextStep / Mpc;
	}
	int getStatus() const {
		return status;
	}
	void setLastStep(double lstep) {
		lastStep = lstep;
	}
	void setNextStep(double nstep) {
		nextStep = nstep;
	}
	void setNextStepMpc(double nstep) {
		nextStep = nstep * Mpc;
	}
	void setStatus(Status stat) {
		status = stat;
	}
};
} // namespace mpc

#endif /* PARTICLE_H_ */

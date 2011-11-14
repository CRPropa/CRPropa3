#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "ThreeVector.h"

// SI units
static const double meter = 1;
static const double second = 1;
static const double kilogram = 1;
static const double ampere = 1;
static const double mol = 1;
static const double kelvin = 1;

// derived units
static const double joule = 1;

// physical constants
static const double eplus = 1.602176487e-19 * ampere * second;
static const double c_light = 2.99792458e+8 * meter / second;
static const double c_squared = c_light * c_light;
static const double amu = 1.660538921e-27 * kilogram;

// other units
static const double electronvolt = eplus * joule;
static const double EeV = 1.e18 * electronvolt;
static const double parsec = 3.0856775807e+16 * meter;
static const double Mpc = 1.e6 * parsec;



class Particle {
private:
	int chargeNumber;
	int massNumber;
	double energy;
	Hep3Vector position;
	Hep3Vector direction;
	double age;
	double step;
	double nextStep;
	int status;

public:
	Particle() {
		chargeNumber = 0;
		massNumber = 0;
		position = Hep3Vector(0);
		direction = Hep3Vector(0);
		energy = 0;
		status = 0;
	}

	// need mass as additional identifier for exitation states in case of nuclei
	void setChargeMassNumber(int Z, int A) {
		chargeNumber = Z;
		massNumber = A;
	}
	void setEnergy(double E) {
		energy = E;
	}
	void setEnergyEeV(double E) {
		energy = E * EeV;
	}
	void setPosition(const Hep3Vector &x) {
		position = x;
	}
	void setPositionMpc(const Hep3Vector &x) {
		position = x * Mpc;
	}
	void setDirection(const Hep3Vector &p) {
		direction = p;
	}
	void setTrajectoryLength(double a) {
		age = a;
	}
	void setTrajectoryLengthMpc(double a) {
		age = a * Mpc;
	}
	void setStep(double s) {
		step = s;
	}
	void setStepMpc(double s) {
		step = s * Mpc;
	}
	void setNextStep(double ns) {
		nextStep = ns;
	}
	void setNextStepMpc(double ns) {
		nextStep = ns * Mpc;
	}
	void setStatus(int s) {
		status = s;
	}

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
		return position * (1 / Mpc);
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

	double getTrajectoryLength() const {
		return age;
	}
	double getTrajectoryLengthMpc() const {
		return age / Mpc;
	}
	double getStep() const {
		return step;
	}
	double getStepMpc() const {
		return step / Mpc;
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
};

#endif /* PARTICLE_H_ */

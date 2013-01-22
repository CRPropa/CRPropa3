#include "mpc/ParticleState.h"

#include <HepPID/ParticleIDMethods.hh>
#include <stdexcept>

namespace mpc {

ParticleState::ParticleState() :
		id(0), pmass(0), energy(0), position(0, 0, 0), direction(-1, 0, 0) {

}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getMag();
}

const Vector3d &ParticleState::getDirection() const {
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
	if (HepPID::isNucleus(id)) {
		pmass = nucleusMass(id);
		charge = HepPID::Z(id) * eplus; // HepPID::charge doesn't work for nuclei
		if (id < 0)
			charge *= -1; // HepPID::Z returns positive charge numbers for anti-nuclei
	} else {
		charge = HepPID::charge(id) * eplus;
	}
}

int ParticleState::getId() const {
	return id;
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getCharge() const {
	return charge;
}

int ParticleState::getChargeNumber() const {
	return HepPID::Z(id);
}

int ParticleState::getMassNumber() const {
	return HepPID::A(id);
}

bool ParticleState::isNucleus() const {
	return HepPID::isNucleus(id);
}

double ParticleState::getLorentzFactor() const {
	if (HepPID::isNucleus(id))
		return energy / (pmass * c_squared);
	else
		throw std::runtime_error(
				"mpc::ParticleState::getLorentzFactor only for nuclei/nucleons");
}

void ParticleState::setLorentzFactor(double gamma) {
	if (HepPID::isNucleus(id))
		energy = gamma * pmass * c_squared;
	else
		throw std::runtime_error(
				"mpc::ParticleState::setLorentzFactor only for nuclei/nucleons");
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3d ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

} // namespace mpc

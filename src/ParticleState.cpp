#include "crpropa/ParticleState.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"

#include "HepPID/ParticleIDMethods.hh"

#include <cstdlib>
#include <cmath>
#include <sstream>

namespace crpropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir): id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.)
{
	setId(id);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir.getUnitVector();
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy);
}

void ParticleState::setId(int newId) {
	id = newId;
	pmass = particleMass(id);
	if (isNucleus(id)) {
		charge = chargeNumber(id) * eplus;
		if (id < 0)
			charge *= -1; // anti-nucleus
	} else {
		charge = HepPID::charge(id) * eplus;
	}
}

double ParticleState::getLorentzFactor() const {
	if (pmass==0)
		return INFINITY;
	return energy/pmass/c_squared + 1;  // should never be inf as long as the energy is representable
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	setEnergy((lf-1) * pmass * c_squared);
}

double ParticleState::getBeta() const {
	return getVelocity().getR2()/c_squared;
}

Vector3d ParticleState::getVelocity() const {
	Vector3d velocity;
	if (pmass==0) 
		velocity = direction*c_light;
	else if (getLorentzFactor()<1.001)  // can happen if if gamma-1 < numericalPrecission
		velocity = direction * sqrt(energy*2/pmass);  // non relativistic case
	else
		velocity = direction * c_light*sqrt(1-1/pow(getLorentzFactor(), 2));

	return velocity;
}

Vector3d ParticleState::getMomentum() const {
	if (pmass==0)
		return direction*energy/c_light;
	else if (getLorentzFactor()<1.001)
		return pmass*getVelocity();
	else
		return getLorentzFactor()*pmass*getVelocity();
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << energy / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "dir = " << direction << ", ";
	ss << "p = " << getMomentum() << " kg*m/s, ";
	ss << "v = " << getVelocity() << " m/s, " ;
	return ss.str();
}

} // namespace crpropa

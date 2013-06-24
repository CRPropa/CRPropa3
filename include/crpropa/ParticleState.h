#ifndef CRPROPA_PARTICLE_STATE_H
#define CRPROPA_PARTICLE_STATE_H

#include "crpropa/Vector3.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"

namespace crpropa {

/**
 @class ParticleState
 @brief State of the particle: ID, energy, position, direction

 The ParticleState defines the state of an ultra-high energy cosmic ray, which
 is assumed to be flying at the exact speed of light. The cosmic ray state is
 uniquely defined by particle ID, energy and position and direction vector.
 For faster lookup mass and charge of the particle are stored as members.
 */
class ParticleState {
private:
	int id; /*< particle ID (Particle Data Group numbering scheme) */
	double energy; /*< total energy */
	Vector3d position; /*< position vector in comoving coordinates */
	Vector3d direction; /*< unit vector of velocity or momentum */

	double pmass;
	double charge;

public:
	ParticleState();
	/* Set position in comoving coordinates */
	void setPosition(const Vector3d &pos);
	/* Get position in comoving coordinates */
	const Vector3d &getPosition() const;

	/* Set direction unit vector. Non unit-vectors are normalized */
	void setDirection(const Vector3d &dir);
	const Vector3d &getDirection() const;

	void setEnergy(double newEnergy);
	double getEnergy() const;

	void setId(int);
	int getId() const;

	/* Electrical charge of the particle */
	double getCharge() const;
	/* Mass of the particle */
	double getMass() const;

	/* Set Lorentz factor and modify the particles energy accordingly */
	void setLorentzFactor(double gamma);
	double getLorentzFactor() const;

	/* Velocity: direction times the speed of light */
	Vector3d getVelocity() const;
	/* Momentum: direction times energy divided by the speed of light */
	Vector3d getMomentum() const;

	/* Check if the particle is a (anti-)nucleus */
	bool isNucleus() const;
	/* If nucleus, return the charge number Z (apositive for anti-nuclei) */
	int getChargeNumber() const;
	/* If nucleus, return the mass number A (also positive for anti-nuclei) */
	int getMassNumber() const;
};

} // namespace crpropa

#endif // CRPROPA_PARTICLE_STATE_H

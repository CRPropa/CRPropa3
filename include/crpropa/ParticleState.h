#ifndef CRPROPA_PARTICLE_STATE_H
#define CRPROPA_PARTICLE_STATE_H

#include "crpropa/Vector3.h"

namespace crpropa {

/**
 @class ParticleState
 @brief State of the particle: ID, energy, position, direction

 The ParticleState defines the state of an ultra-high energy cosmic ray, which
 is assumed to be traveling at the exact speed of light.
 The cosmic ray state is defined by particle ID, energy and position and
 direction vector.
 For faster lookup mass and charge of the particle are stored as members.
 */
class ParticleState {
private:
	int id; ///< particle ID (Particle Data Group numbering scheme)
	double energy; ///< total energy
	Vector3d position; ///< position vector in comoving coordinates
	Vector3d direction; ///< unit vector of velocity or momentum
	double pmass; ///< particle rest mass
	double charge; ///< particle charge

public:
	ParticleState(int id = 0, double energy = 0,
			Vector3d position = Vector3d(0, 0, 0),
			Vector3d direction = Vector3d(-1, 0, 0));

	/// Set position in comoving coordinates
	void setPosition(const Vector3d &pos);
	/// Get position in comoving coordinates
	const Vector3d &getPosition() const;

	/// Set direction unit vector, non unit-vectors are normalized
	void setDirection(const Vector3d &dir);
	/// Get direction unit vector
	const Vector3d &getDirection() const;

	/// Set energy in [J]
	void setEnergy(double newEnergy);
	/// Get energy in [J]
	double getEnergy() const;

	/// Set particle ID
	void setId(int);
	/// Get particle ID
	int getId() const;

	std::string getDescription() const;

	// ======== Helper methods ========

	/// Electrical charge of the particle in [A]
	double getCharge() const;
	/// Mass of the particle in [kg]
	double getMass() const;

	/// Set Lorentz factor and modify the particle's energy accordingly
	void setLorentzFactor(double gamma);
	/// Get Lorentz factor
	double getLorentzFactor() const;

	/// Velocity: direction times the speed of light in [m/s]
	Vector3d getVelocity() const;
	/// Momentum: direction times energy divided by the speed of light [kg m/s]
	Vector3d getMomentum() const;
};

} // namespace crpropa

#endif // CRPROPA_PARTICLE_STATE_H

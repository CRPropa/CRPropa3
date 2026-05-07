#ifndef CRPROPA_PARTICLE_STATE_H
#define CRPROPA_PARTICLE_STATE_H

#include "crpropa/Vector3.h"

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

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
	/** Constructor for a particle state.
	 @param id			id of the particle following the PDG numbering scheme
	 @param energy		energy of the particle [in Joules]
	 @param position	vector containing the coordinates of the particle [in meters]
	 @param direction	vector containing the direction of motion of the particle
	 */
	ParticleState(int id = 0, double energy = 0,
			Vector3d position = Vector3d(0, 0, 0),
			Vector3d direction = Vector3d(-1, 0, 0));

	/** Set particle position.
	 In simulations including cosmological effects, the position is given in comoving coordinates.
	 @param pos		vector containing the coordinates of the particle [in meters]
	*/
	void setPosition(const Vector3d &pos);
	/** Get position of particle.
	 @returns Position vector of particle. If cosmological effects are included, the coordinates are comoving.
	 */
	const Vector3d &getPosition() const;

	/** Set direction unit vector, non unit-vectors are normalized
	 @param dir	vector containing the direction of motion of the particle
	 */
	void setDirection(const Vector3d &dir);
	/** Get direction unit vector
	 @returns Normalized vector containing direction of motion of particle.
	 */
	const Vector3d &getDirection() const;

	/** Set energy of particle.
	 @param newEnergy	energy to be assigned to particle [in Joules]
	 */
	void setEnergy(double newEnergy);
	/** Get energy of particle.
	 @returns Energy of particle [in Joules]
	 */
	double getEnergy() const;
	/** Get rigidity of particle, defined as E/(Z*e).
	 @returns Rigidity of the particle [in Volts]
	 */
	double getRigidity() const;

	/** Set particle ID.
	 This follows the PDG numbering scheme:
	  https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
	 @param newId		id to be assigned to the particle 
	 */
	void setId(int newId);
	/** Get particle ID
	 @returns Particle ID (in PDG format).
	 */
	int getId() const;

	std::string getDescription() const;

	// ======== Helper methods ========

	/** Get electrical charge of the particle.
	 @returns Charge of the particle [in Coulombs]
	 */
	double getCharge() const;
	/** Get mass of the particle.
	 @returns Mass of the particle [kg]
	 */
	double getMass() const;

	/** Set Lorentz factor and modify the particle's energy accordingly.
	 @param gamma		Lorentz factor
	 */
	void setLorentzFactor(double gamma);
	/** Get Lorentz factor
	 @returns Lorentz factor of particle
	 */
	double getLorentzFactor() const;

	/** Get velocity: direction times the speed of light.
	 @returns Velocity of particle [m/s]
	 */
	Vector3d getVelocity() const;
	/** Get momentum: direction times energy divided by the speed of light 
	 @returns The momentum [kg m/s]
	*/
	Vector3d getMomentum() const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PARTICLE_STATE_H

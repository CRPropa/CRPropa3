#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"
#include "mpc/Common.h"
#include "mpc/Nucleus.h"

namespace mpc {

/**
 @class ParticleState
 @brief State of the particle: ID, energy, position, direction
 */
class ParticleState {
private:
	int id;
	double pmass;
	double energy;
	Vector3d position;
	Vector3d direction;

public:
	ParticleState();
	void setPosition(const Vector3d &pos);
	const Vector3d &getPosition() const;

	void setDirection(const Vector3d &dir);
	const Vector3d &getDirection() const;

	void setEnergy(const double newEnergy);
	double getEnergy() const;

	void setId(const int);
	int getId() const;

	int getChargeNumber() const;
	double getCharge() const;

	int getMassNumber() const;
	double getMass() const;

	// convenience functions
	bool isNucleus() const;

	void setLorentzFactor(const double gamma);
	double getLorentzFactor() const;

	Vector3d getVelocity() const;
	Vector3d getMomentum() const;
};

} // namespace mpc

#endif /* PARTICLE_H_ */

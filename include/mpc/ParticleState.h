#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"
#include "mpc/HepPID.h"

namespace mpc {

/**
 @class ParticleState
 @brief The state defining the particle.
 */
class ParticleState {
public:
	void setPosition(const Vector3 &pos);
	const Vector3 &getPosition() const;

	void setDirection(const Vector3 &dir);
	const Vector3 &getDirection() const;

	void setEnergy(double newEnergy);
	double getEnergy() const;

	void setId(int);
	int getId() const;

	int getChargeNumber() const;
	double getCharge() const;

	int getMassNumber() const;
	double getMass() const;

	// convenience functions
	double getLorentzFactor() const;
	Vector3 getVelocity() const;
	Vector3 getMomentum() const;

private:
	int id;
	double energy;
	Vector3 position;
	Vector3 direction;
};

} // namespace mpc

#endif /* PARTICLE_H_ */

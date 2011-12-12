#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"

namespace mpc {

class ParticleState {
public:
	const Vector3 &getPosition() const;
	void setPosition(const Vector3 &pos);

	const Vector3 &getDirection() const;
	void setDirection(const Vector3 &dir);

	int getId() const;
	void setId(int);

	double getEnergy() const;
	void setEnergy(double newEnergy);

	double getChargeNumber() const;

	double getMassNumber() const;

	double getMass() const;

	// convenience functions
	double getLorentzFactor() const;
	Vector3 getVelocity() const;
	Vector3 getMomentum() const;

private:
	double energy;
	Vector3 position;
	Vector3 direction;
	double mass;
	int id;

	void setMass(double newMass);
};

} // namespace mpc

#endif /* PARTICLE_H_ */

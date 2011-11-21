#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "mpc/Vector3.h"
#include "mpc/Units.h"

namespace mpc {

class ParticleState {
public:
	enum Type {
		Gamma, Lepton, Hadron
	};

	union {
		size_t baryonNumber;
		size_t leptonNumber;
		size_t massNumber;
	};

	double getEnergy() const;
	void setEnergy(double newEnergy);

	const Vector3 &getPosition() const;
	void setPosition(const Vector3 &pos);

	const Vector3 &getDirection() const;
	void setDirection(const Vector3 &dir);

	double getChargeNumber() const;
	void setChargeNumber(size_t charge);

	double getMass() const;
	void setMass(double newMass);

	Type getType();
	void setType(Type t);

	// convenience
	double getLorentzFactor() const;
	Vector3 getVelocity() const;
	Vector3 getMomentum() const;

private:
	double energy;
	Vector3 position;
	Vector3 direction;
	size_t chargeNumber;
	double mass;
	Type type;

};

} // namespace mpc

#endif /* PARTICLE_H_ */

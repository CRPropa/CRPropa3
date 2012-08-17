#ifndef MPC_OBSERVER_H_
#define MPC_OBSERVER_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class SmallObserverSphere
 @brief Flags particles when entering the sphere.

 Particles are detected when the current position is inside and the previous position outside the sphere.
 In this case the candidate is by default flagged "Detected" made inactive.
 This module limits the next step size to prevent candidates from overshooting.
 */
class SmallObserverSphere: public Module {
private:
	Vector3d center;
	double radius;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	void updateDescription();

public:
	SmallObserverSphere(Vector3d center, double radius, std::string flag =
			"Detected", std::string flagValue = "", bool makeInactive = true);
	void process(Candidate *candidate) const;
};

/**
 @class LargeObserverSphere
 @brief Flags particles when leaving the sphere.

 Particles are detected when the current position is outside and the previous position inside the sphere.
 In this case the candidate is by default flagged "Detected" made inactive.
 This module limits the next step size to prevent candidates from overshooting.
 */
class LargeObserverSphere: public Module {
private:
	Vector3d center;
	double radius;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	void updateDescription();

public:
	LargeObserverSphere(Vector3d center, double radius, std::string flag =
			"Detected", std::string flagValue = "", bool makeInactive = true);
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* MPC_OBSERVER_H_ */

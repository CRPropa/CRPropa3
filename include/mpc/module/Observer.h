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

public:
	SmallObserverSphere(Vector3d center = Vector3d(0.), double radius = 0,
			std::string flag = "Detected", std::string flagValue = "",
			bool makeInactive = true);
	void process(Candidate *candidate) const;
	void setCenter(Vector3d center);
	void setRadius(double radius);
	void setFlag(std::string flag, std::string flagValue);
	void setMakeInactive(bool makeInactive);
	std::string getDescription() const;
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
	LargeObserverSphere(Vector3d center = Vector3d(0.), double radius = 0,
			std::string flag = "Detected", std::string flagValue = "",
			bool makeInactive = true);
	void process(Candidate *candidate) const;
	void setCenter(Vector3d center);
	void setRadius(double radius);
	void setFlag(std::string flag, std::string flagValue);
	void setMakeInactive(bool makeInactive);
	std::string getDescription() const;
};

/**
 @class Observer1D
 @brief Observer for 1D simulations

 Particles are detected once their x-position gets smaller than a 0.
 In this case the candidate is by flagged "Detected" made inactive.
 This module limits the next step size to prevent candidates from overshooting.
 */
class Observer1D: public Module {
public:
	Observer1D();
	void process(Candidate *candidate) const;
};

/**
 @class DetectAll
 @brief Flags all

 Module that flags every particle it processes.
 Optionally inactivates the particle.
 */
class DetectAll: public Module {
private:
	std::string flag;
	std::string flagValue;
	bool makeInactive;

public:
	DetectAll(std::string flag = "Detected", std::string flagValue = "",
			bool makeInactive = true);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_OBSERVER_H_ */

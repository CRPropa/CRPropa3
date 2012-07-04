#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_

#include "mpc/Module.h"
#include "mpc/Vector3.h"
#include <sstream>

namespace mpc {

/**
 @class MaximumTrajectoryLength
 @brief Stops propagation after a maximum trajectory length.

 This modules deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is no exceeded.
 */
class MaximumTrajectoryLength: public Module {
private:
	double maxLength;

public:
	MaximumTrajectoryLength(double length = 0);
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergy
 @brief Stops propagation if energy drops under the set minimum energy.

 This modules deactivates the candidate below a give minimum energy.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumEnergy: public Module {
private:
	double minEnergy;

public:
	MinimumEnergy(double minEnergy = 0);
	void setMinimumEnergy(double energy);
	double getMinimumEnergy() const;
	void process(Candidate *candidate) const;
};

/**
 @class SmallObserverSphere
 @brief Detects particles when entering the sphere.
 */
class SmallObserverSphere: public Module {
	void updateDescription();
public:
	double radius;
	Vector3d center;
	std::string flag;
	std::string flagValue;
	bool makeInactive;

	SmallObserverSphere(Vector3d center, double radius, std::string flag =
			"Detected", std::string flagValue = "");
	void setMakeInactive(bool makeInactive);
	void process(Candidate *candidate) const;
};

/**
 @class PeriodicBox
 @brief Box with periodic boundaries

 If a particle leaves the periodic box it is placed at the opposite side and its initial (source) position changed accordingly.
 Particles can overshoot (be outside of the box during the step) since the step size is not limited by this module.
 */
class PeriodicBox: public Module {
public:
	Vector3d origin;
	Vector3d size;

	PeriodicBox(Vector3d origin, Vector3d size);
	void process(Candidate *candidate) const;
};

/**
 @class CubicBoundary
 @brief Flags a particle when exiting the cube.
 */
class CubicBoundary: public Module {
protected:
	Vector3d origin;
	double size;
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;
	void updateDescription();
public:
	CubicBoundary(Vector3d origin, double size, std::string flag = "OutOfBounds",
			std::string flagValue = "");
	void setMakeInactive(bool makeInactive);
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

/**
 @class SphericalBoundary
 @brief Flag a particle when leaving the sphere.
 */
class SphericalBoundary: public Module {
protected:
	Vector3d center;
	double radius;
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;
	void updateDescription();
public:
	SphericalBoundary(Vector3d center, double radius, std::string flag =
			"OutOfBounds", std::string flagValue = "");
	void setMakeInactive(bool makeInactive);
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

/**
 @class EllipsoidalBoundary
 @brief Flags a particle when leaving the ellipsoid.
 */
class EllipsoidalBoundary: public Module {
protected:
	Vector3d focalPoint1;
	Vector3d focalPoint2;
	double majorAxis;
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;
	void updateDescription();
public:
	EllipsoidalBoundary(Vector3d focalPoint1, Vector3d focalPoint2,
			double majorAxis, std::string flag = "OutOfBounds",
			std::string flagValue = "");
	void setMakeInactive(bool makeInactive);
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

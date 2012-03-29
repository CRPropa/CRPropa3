#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_

#include "mpc/Module.h"
#include "mpc/Vector3.h"
#include <sstream>

namespace mpc {

/**
 @class MaximumTrajectoryLength
 @brief Stops propagation after the set maximum trajectory length.
 */
class MaximumTrajectoryLength: public Module {
private:
	double maxLength;

public:
	MaximumTrajectoryLength(double length);
	void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergy
 @brief Stops propagation if energy drops under the set minimum energy.
 */
class MinimumEnergy: public Module {
public:
	double minEnergy;

	MinimumEnergy(double minEnergy);
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
	Vector3 center;
	std::string flag;
	std::string flagValue;
	bool makeInactive;

	SmallObserverSphere(Vector3 center, double radius, std::string flag =
			"Detected", std::string flagValue = "");
	void setMakeInactive(bool makeInactive);
	void process(Candidate *candidate) const;
};

/**
 @class CubicBoundary
 @brief Flags a particle when exiting the cube.
 */
class CubicBoundary: public Module {
protected:
	Vector3 origin;
	double size;
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;
	void updateDescription();
public:
	CubicBoundary(Vector3 origin, double size, std::string flag = "OutOfBounds",
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
	Vector3 center;
	double radius;
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;
	void updateDescription();
public:
	SphericalBoundary(Vector3 center, double radius, std::string flag =
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
	Vector3 focalPoint1;
	Vector3 focalPoint2;
	double majorAxis;
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;
	void updateDescription();
public:
	EllipsoidalBoundary(Vector3 focalPoint1, Vector3 focalPoint2,
			double majorAxis, std::string flag = "OutOfBounds",
			std::string flagValue = "");
	void setMakeInactive(bool makeInactive);
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

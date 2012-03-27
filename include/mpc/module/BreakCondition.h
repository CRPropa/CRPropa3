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
	MaximumTrajectoryLength(double length) :
			maxLength(length) {
	}

	void process(Candidate *candidate) const {
		double l = candidate->getTrajectoryLength();
		if (l >= maxLength)
			candidate->setActive(false);
		else
			candidate->limitNextStep(maxLength - l);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc";
		return s.str();
	}
};

/**
 @class MinimumEnergy
 @brief Stops propagation if energy drops under the set minimum energy.
 */
class MinimumEnergy: public Module {
public:
	double minEnergy;

	MinimumEnergy(double minEnergy) :
			minEnergy(minEnergy) {
	}

	void process(Candidate *candidate) const {
		if (candidate->current.getEnergy() <= minEnergy)
			candidate->setActive(false);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Minimum energy: " << minEnergy / EeV << " EeV";
		return s.str();
	}
};

/**
 @class SmallObserverSphere
 @brief Detects particles when entering the sphere.
 */
class SmallObserverSphere: public Module {
public:
	double radius;
	Vector3 center;

	SmallObserverSphere(Vector3 center, double radius) :
			center(center), radius(radius) {
	}

	void process(Candidate *candidate) const {
		double d = (candidate->current.getPosition() - center).mag();
		if (d <= radius * 1.01) {
			candidate->setActive(false);
			candidate->setProperty("Detected", "");
		} else
			candidate->limitNextStep((d - radius));
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Small observer sphere: " << radius << " radius around " << center;
		return s.str();
	}
};

/**
 @class CubicBoundary
 @brief Flags a particle when exiting the cube.
 */
class CubicBoundary: public Module {
protected:
	Vector3 origin;
	double size;
	bool limitStep;
	double margin;

public:
	CubicBoundary(Vector3 origin, double size) :
			origin(origin), size(size), limitStep(false) {
	}

	CubicBoundary(Vector3 origin, double size, double margin) :
			origin(origin), size(size), limitStep(true), margin(margin) {
	}

	void process(Candidate *candidate) const {
		Vector3 relPos = candidate->current.getPosition() - origin;
		double lo = std::min(relPos.x(), std::min(relPos.y(), relPos.z()));
		double hi = std::max(relPos.x(), std::max(relPos.y(), relPos.z()));
		if ((lo <= 0.) or (hi >= size))
			candidate->setActive(false);
		if (limitStep) {
			candidate->limitNextStep(lo + margin);
			candidate->limitNextStep(size - hi + margin);
		}
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Cubic Boundary: origin " << origin << ", size " << size;
		return s.str();
	}
};

/**
 @class SphericalBoundary
 @brief Flags a particle when exiting the sphere.
 */
class SphericalBoundary: public Module {
protected:
	Vector3 center;
	double radius;
	bool limitStep;
	double margin;

public:
	SphericalBoundary(Vector3 center, double radius) :
			center(center), radius(radius), limitStep(false) {
	}

	SphericalBoundary(Vector3 center, double radius, double margin) :
			center(center), radius(radius), limitStep(true), margin(margin) {
	}

	void process(Candidate *candidate) const {
		double d = (candidate->current.getPosition() - center).mag();
		if (d >= radius)
			candidate->setActive(false);
		if (limitStep)
			candidate->limitNextStep(radius - d + margin);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Spherical Boundary: radius " << radius << " around " << center;
		return s.str();
	}
};

/**
 @class EllipsoidalBoundary
 @brief Flags a particle when leaving then ellipsoid.
 */
class EllipsoidalBoundary: public Module {
protected:
	Vector3 focalPoint1;
	Vector3 focalPoint2;
	double majorAxis;
	bool limitStep;
	double margin;

public:
	EllipsoidalBoundary(Vector3 focalPoint1, Vector3 focalPoint2,
			double majorAxis) {
		this->focalPoint1 = focalPoint1;
		this->focalPoint2 = focalPoint2;
		this->majorAxis = majorAxis;
		this->limitStep = false;
	}

	EllipsoidalBoundary(Vector3 focalPoint1, Vector3 focalPoint2,
			double majorAxis, double margin) {
		this->focalPoint1 = focalPoint1;
		this->focalPoint2 = focalPoint2;
		this->majorAxis = majorAxis;
		this->limitStep = true;
		this->margin = margin;
	}

	void process(Candidate *candidate) const {
		Vector3 pos = candidate->current.getPosition();
		double d = (pos - focalPoint1).mag() + (pos - focalPoint2).mag();
		if (d >= majorAxis)
			candidate->setActive(false);
		if (limitStep)
			candidate->limitNextStep(majorAxis - d + margin);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Ellipsoidal Boundary: F1 = " << focalPoint1 << ", F2 = "
				<< focalPoint2 << ", major axis = " << majorAxis;
		return s.str();
	}
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

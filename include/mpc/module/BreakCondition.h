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
	std::string flag;
	std::string flagValue;
	bool makeInactive;

	SmallObserverSphere(Vector3 center, double radius, std::string flag =
			"Detected", std::string flagValue = "") {
		this->center = center;
		this->radius = radius;
		this->flag = flag;
		this->flagValue = flagValue;
		this->makeInactive = true;
	}

	void setMakeInactive(bool makeInactive) {
		this->makeInactive = makeInactive;
	}

	void process(Candidate *candidate) const {
		double d = (candidate->current.getPosition() - center).mag();
		if (d <= radius * 1.01) {
			candidate->addProperty("Detected", "");
			if (makeInactive)
				candidate->setActive(false);
		}
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
	double margin;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	bool limitStep;

public:
	CubicBoundary(Vector3 origin, double size, std::string flag = "OutOfBounds",
			std::string flagValue = "") {
		this->origin = origin;
		this->size = size;
		this->flag = flag;
		this->flagValue = flagValue;
		this->makeInactive = false;
		this->limitStep = false;
		this->margin = 0;
	}

	void setMakeInactive(bool makeInactive) {
		this->makeInactive = makeInactive;
	}

	void setLimitStep(bool limitStep, double margin) {
		this->limitStep = limitStep;
		this->margin = margin;
	}

	void process(Candidate *candidate) const {
		Vector3 relPos = candidate->current.getPosition() - origin;
		double lo = std::min(relPos.x(), std::min(relPos.y(), relPos.z()));
		double hi = std::max(relPos.x(), std::max(relPos.y(), relPos.z()));
		if ((lo <= 0.) or (hi >= size)) {
			candidate->addProperty(flag, flagValue);
			if (makeInactive)
				candidate->setActive(false);
		}
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

public:
	SphericalBoundary(Vector3 center, double radius, std::string flag =
			"OutOfBounds", std::string flagValue = "") {
		this->center = center;
		this->radius = radius;
		this->flag = flag;
		this->flagValue = flagValue;
		this->makeInactive = false;
		this->limitStep = false;
		this->margin = 0;
	}

	void setMakeInactive(bool makeInactive) {
		this->makeInactive = makeInactive;
	}

	void setLimitStep(bool limitStep, double margin) {
		this->limitStep = limitStep;
		this->margin = margin;
	}

	void process(Candidate *candidate) const {
		double d = (candidate->current.getPosition() - center).mag();
		if (d >= radius) {
			candidate->addProperty(flag, flagValue);
			if (makeInactive)
				candidate->setActive(false);
		}
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

public:
	EllipsoidalBoundary(Vector3 focalPoint1, Vector3 focalPoint2,
			double majorAxis, std::string flag = "OutOfBounds",
			std::string flagValue = "") {
		this->focalPoint1 = focalPoint1;
		this->focalPoint2 = focalPoint2;
		this->majorAxis = majorAxis;
		this->flag = flag;
		this->flagValue = flagValue;
		this->makeInactive = false;
		this->limitStep = false;
		this->margin = 0;
	}

	void setMakeInactive(bool makeInactive) {
		this->makeInactive = makeInactive;
	}

	void setLimitStep(bool limitStep, double margin) {
		this->limitStep = limitStep;
		this->margin = margin;
	}

	void process(Candidate *candidate) const {
		Vector3 pos = candidate->current.getPosition();
		double d = (pos - focalPoint1).mag() + (pos - focalPoint2).mag();
		if (d >= majorAxis) {
			candidate->addProperty(flag, flagValue);
			if (makeInactive)
				candidate->setActive(false);
		}
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

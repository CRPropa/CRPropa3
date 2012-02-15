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
public:
	double maxLength;

	MaximumTrajectoryLength(double maxLength) {
		this->maxLength = maxLength;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		if (candidate->getTrajectoryLength() >= maxLength)
			candidate->setStatus(Candidate::Stopped);
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

	MinimumEnergy(double minEnergy) {
		this->minEnergy = minEnergy;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		if (candidate->current.getEnergy() <= minEnergy)
			candidate->setStatus(Candidate::Stopped);
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

	SmallObserverSphere(Vector3 center, double radius) {
		this->radius = radius;
		this->center = center;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		double d = (candidate->current.getPosition() - center).mag();
		if (d <= radius)
			candidate->setStatus(Candidate::Detected);
		else
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
	Candidate::Status flag;
	bool hardBoundary;
	double margin;

public:
	CubicBoundary(Vector3 origin, double size, Candidate::Status flag) {
		this->origin = origin;
		this->size = size;
		this->flag = flag;
	}

	CubicBoundary(Vector3 origin, double size, double margin,
			Candidate::Status flag) {
		this->origin = origin;
		this->size = size;
		this->flag = flag;
		this->hardBoundary = true;
		this->margin = margin;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		Vector3 relPos = candidate->current.getPosition() - origin;
		double lo = std::min(relPos.x(), std::min(relPos.y(), relPos.z()));
		double hi = std::max(relPos.x(), std::max(relPos.y(), relPos.z()));
		if ((lo <= 0.) && (hi >= size))
			candidate->setStatus(flag);
		if (hardBoundary) {
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
	Candidate::Status flag;
	bool hardBoundary;
	double margin;

public:
	SphericalBoundary(Vector3 center, double radius, Candidate::Status flag) {
		this->center = center;
		this->radius = radius;
		this->flag = flag;
	}

	SphericalBoundary(Vector3 center, double radius, double margin,
			Candidate::Status flag) {
		this->center = center;
		this->radius = radius;
		this->flag = flag;
		this->hardBoundary = true;
		this->margin = margin;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		double d = (candidate->current.getPosition() - center).mag();
		if (d >= radius)
			candidate->setStatus(flag);
		if (hardBoundary)
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
	Candidate::Status flag;
	bool hardBoundary;
	double margin;

public:
	EllipsoidalBoundary(Vector3 focalPoint1, Vector3 focalPoint2,
			double majorAxis, Candidate::Status flag) {
		this->focalPoint1 = focalPoint1;
		this->focalPoint2 = focalPoint2;
		this->majorAxis = majorAxis;
		this->flag = flag;
	}

	EllipsoidalBoundary(Vector3 focalPoint1, Vector3 focalPoint2,
			double majorAxis, double margin, Candidate::Status flag) {
		this->focalPoint1 = focalPoint1;
		this->focalPoint2 = focalPoint2;
		this->majorAxis = majorAxis;
		this->flag = flag;
		this->hardBoundary = true;
		this->margin = margin;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		Vector3 pos = candidate->current.getPosition();
		double d = (pos - focalPoint1).mag() + (pos - focalPoint2).mag();
		if (d >= majorAxis)
			candidate->setStatus(flag);
		if (hardBoundary)
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

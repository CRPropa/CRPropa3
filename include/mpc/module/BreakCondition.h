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

	void process(Candidate *candidate) {
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

	void process(Candidate *candidate) {
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
 @class LargeObserverSphere
 @brief Detects particles when leaving the sphere.
 */
class LargeObserverSphere: public Module {
public:
	double radius;
	Vector3 center;

	LargeObserverSphere(Vector3 center, double radius) {
		this->radius = radius;
		this->center = center;
	}

	void process(Candidate *candidate) {
		double d = (candidate->current.getPosition() - center).mag();
		if (d >= radius)
			candidate->setStatus(Candidate::Detected);
		else
			candidate->limitNextStep(radius - d);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Large observer sphere: " << radius / Mpc << " Mpc radius around "
				<< center / Mpc << " Mpc";
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

	void process(Candidate *candidate) {
		double d = (candidate->current.getPosition() - center).mag();
		if (d <= radius)
			candidate->setStatus(Candidate::Detected);
		else
			candidate->limitNextStep((d - radius));
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Small observer sphere: " << radius / Mpc << " Mpc radius around "
				<< center / Mpc << " Mpc";
		return s.str();
	}
};

/**
 @class SimulationBox
 @brief Declares particle out of bounds when exiting the simulation box.
 */
class SimulationBox: public Module {
public:
	Vector3 origin;
	double size;
	double margin;

	SimulationBox(Vector3 origin, double size, double margin) {
		this->origin = origin;
		this->size = size;
		this->margin = margin;
	}

	void process(Candidate *candidate) {
		Vector3 relPos = candidate->current.getPosition() - origin;
		double lo = std::min(relPos.x(), std::min(relPos.y(), relPos.z()));
		double hi = std::max(relPos.x(), std::max(relPos.y(), relPos.z()));
		if (lo < 0.)
			candidate->setStatus(Candidate::OutOfBounds);
		if (hi > size)
			candidate->setStatus(Candidate::OutOfBounds);
		candidate->limitNextStep(lo + margin);
		candidate->limitNextStep(size - hi + margin);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Simulation Box: " << origin / Mpc << " - "
				<< origin + Vector3(size, size, size) << " Mpc, margin "
				<< margin << " Mpc";
		return s.str();
	}
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

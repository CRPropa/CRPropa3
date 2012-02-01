#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_

#include "mpc/Module.h"
#include "mpc/Vector3.h"
#include <sstream>

namespace mpc {

class MaximumTrajectoryLength: public Module {
public:
	double maxLength;

	MaximumTrajectoryLength(double maxLength) {
		this->maxLength = maxLength;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		if (candidate->getTrajectoryLength() >= maxLength)
			candidate->setStatus(Candidate::ReachedMaxTrajectoryLength);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc";
		return s.str();
	}
};

class MinimumEnergy: public Module {
public:
	double minEnergy;

	MinimumEnergy(double minEnergy) {
		this->minEnergy = minEnergy;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		if (candidate->current.getEnergy() <= minEnergy)
			candidate->setStatus(Candidate::BelowEnergyThreshold);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Minimum energy: " << minEnergy / EeV << " EeV";
		return s.str();
	}
};

class LargeObserverSphere: public Module {
public:
	double radius;
	Vector3 center;

	LargeObserverSphere(double radius, Vector3 center) {
		this->radius = radius;
		this->center = center;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		double d = (candidate->current.getPosition() - center).mag();
		if (d >= radius)
			candidate->setStatus(Candidate::Detected);
		candidate->limitNextStep(radius - d);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Large observer sphere: " << radius / Mpc << " Mpc radius around "
				<< center / Mpc << " Mpc";
		return s.str();
	}
};

class SmallObserverSphere: public Module {
public:
	double radius;
	Vector3 center;

	SmallObserverSphere(double radius, Vector3 center) {
		this->radius = radius;
		this->center = center;
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		double d = (candidate->current.getPosition() - center).mag();
		if (d <= radius)
			candidate->setStatus(Candidate::Detected);
		candidate->limitNextStep((d - radius));
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "Small observer sphere: " << radius / Mpc << " Mpc radius around "
				<< center / Mpc << " Mpc";
		return s.str();
	}
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

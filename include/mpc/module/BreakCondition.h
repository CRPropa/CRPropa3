#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_

#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/Module.h"

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
			candidate->setStatus(Candidate::ReachedMaxTime);
	}

	std::string getDescription() const {
		std::stringstream s;
		s << "MaximumTrajectoryLength: " << maxLength / Mpc << " Mpc";
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
		s << "MinimumEnergy: " << minEnergy / EeV << " EeV";
		return s.str();
	}
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_

#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/Propagator.h"

#include <sstream>

namespace mpc {

class MaximumTrajectoryLength: public Feature {
public:
	double maxLength;

	MaximumTrajectoryLength(double maxLength) {
		this->maxLength = maxLength;
	}

	void apply(Candidate &candidate, size_t priority) {
		if (candidate.getTrajectoryLength() >= maxLength)
			candidate.setStatus(Candidate::ReachedMaxTime);
	}

	std::string description() const {
		std::stringstream s;
		s << "MaximumTrajectoryLength " << maxLength;
		return s.str();
	}
};

class MinimumEnergy: public Feature {
public:
	double minEnergy;

	MinimumEnergy(double minEnergy) {
		this->minEnergy = minEnergy;
	}

	void apply(Candidate &candidate, size_t priority) {
		if (candidate.next.getEnergy() <= minEnergy)
			candidate.setStatus(Candidate::BelowEnergyThreshold);
	}

	std::string description() const {
		std::stringstream s;
		s << "MinimumEnergy " << minEnergy;
		return s.str();
	}
};

} // namespace mpc

#endif /* BREAKCONDITION_H_ */

/*
 * BreakCondition.h
 *
 *  Created on: Nov 11, 2011
 *      Author: walz
 */

#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_
namespace mpc {
class MaximumTrajectoryLength {
public:
	double maxLength;

	MaximumTrajectoryLength(double maxLength) {
		this->maxLength = maxLength;
	}

	void apply(Candidate &candidate) {
		if (candidate.getTrajectoryLength() >= maxLength)
			candidate.setStatus(Candidate::AboveMaxTime);
	}
};

class MinimumEnergy {
public:
	double minEnergy;

	MinimumEnergy(double minEnergy) {
		this->minEnergy = minEnergy;
	}

	void apply(Candidate &candidate) {
		if (candidate.current.getEnergy() <= minEnergy)
			candidate.setStatus(Candidate::BelowEnergyThreshold);
	}
};
} // namespace mpc

#endif /* BREAKCONDITION_H_ */

/*
 * BreakCondition.h
 *
 *  Created on: Nov 11, 2011
 *      Author: walz
 */

#ifndef BREAKCONDITION_H_
#define BREAKCONDITION_H_

class MaximumTrajectoryLength {
public:
	double maxLength;

	MaximumTrajectoryLength(double maxLength) {
		this->maxLength = maxLength;
	}

	void apply(Particle &particle) {
		if (particle.getTrajectoryLength() >= maxLength)
			particle.setStatus(1);
	}
};

class MinimumEnergy {
public:
	double minEnergy;

	MinimumEnergy(double minEnergy) {
		this->minEnergy = minEnergy;
	}

	void apply(Particle &particle) {
		if (particle.getEnergy() <= minEnergy)
			particle.setStatus(1);
	}
};


#endif /* BREAKCONDITION_H_ */

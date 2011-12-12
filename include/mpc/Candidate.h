#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include "mpc/ParticleState.h"

#include <vector>

namespace mpc {

class Candidate {
public:
	enum Status {
		Active = 0,
		Detected,
		ReachedMaxTime,
		BelowEnergyThreshold,
		Decayed,
		ObserverNotReachable,
		UserDefined
	};

	ParticleState current;
	ParticleState last;
	ParticleState initial;
	std::vector< Candidate > secondaries;

	Candidate();

	double getTrajectoryLength() const;
	void setTrajectoryLength(double a);

	double getCurrentStep() const;
	void setCurrentStep(double lstep);

	double getNextStep() const;
	void setNextStep(double nstep);

	void setNextStepUpperLimit(double upper);
	void setNextStepLowerLimit(double lower);

	Status getStatus() const;
	void setStatus(Status stat);
private:
	double age;
	double currentStep, nextStep; // stepsize in [m]
	Status status;
};

} // namespace mpc

#endif /* CANDIDATE_H_ */

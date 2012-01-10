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
		ReachedMaxTrajectoryLength,
		BelowEnergyThreshold,
		ObserverNotReachable,
		UserDefined
	};

	ParticleState current;
	ParticleState last;
	ParticleState initial;
	std::vector< Candidate > secondaries;

	Candidate();

	double getRedshift() const;
	void setRedshift(double z);

	double getTrajectoryLength() const;
	void setTrajectoryLength(double a);

	double getCurrentStep() const;
	void setCurrentStep(double lstep);

	double getNextStep() const;
	void setNextStep(double step);
	void limitNextStep(double step);

	Status getStatus() const;
	void setStatus(Status stat);

private:
	double redshift, trajectoryLength;
	double currentStep, nextStep; // stepsize in [m]
	Status status;
};

} // namespace mpc

#endif /* CANDIDATE_H_ */

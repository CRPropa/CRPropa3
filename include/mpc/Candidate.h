#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include "mpc/ParticleState.h"
#include "mpc/SharedPointer.h"

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
	std::vector< shared_ptr<Candidate> > secondaries;

	Candidate();

	double getTrajectoryLength() const;
	void setTrajectoryLength(double a);

	double getCurrentStep() const;
	void setCurrentStep(double lstep);

	double getNextStep() const;
	void setNextStep(double nstep);

	Status getStatus() const;
	void setStatus(Status stat);
private:
	double age;
	double currentStep, nextStep;
	Status status;
};

} // namespace mpc

#endif /* CANDIDATE_H_ */

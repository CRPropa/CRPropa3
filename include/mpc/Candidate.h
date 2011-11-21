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

	ParticleState next;
	ParticleState last;
	ParticleState initial;
	std::vector< shared_ptr<Candidate> > secondaries;

	Candidate();

	double getTrajectoryLength() const;
	void setTrajectoryLength(double a);

	double getLastStep() const;
	void setLastStepSize(double lstep);

	double getNextStep() const;
	void setNextStepSize(double nstep);

	Status getStatus() const;
	void setStatus(Status stat);
private:
	double age;
	double lastStepSize, nextStepSize;
	Status status;
};

} // namespace mpc

#endif /* CANDIDATE_H_ */

#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include "mpc/Particle.h"

namespace mpc {

class Candidate {
public:
	enum Status {
		Active = 0,
		Detected,
		ReachedMaxTime,
		BelowEnergyThreshold,
		Decayed,
		ObserverNotReachable
	};

	Particle current;
	Particle initial;
	Candidate *parent;

	Candidate();

	double getTrajectoryLength() const;

	double getLastStep() const;
	void setLastStep(double lstep);

	double getNextStep() const;
	void setNextStep(double nstep);

	Status getStatus() const;
	void setStatus(Status stat);
private:
	double age;
	double lastStep, nextStep;
	Status status;
};

Candidate::Candidate() :
		parent(0), age(0), lastStep(0), nextStep(0), status(Active) {

}

double Candidate::getTrajectoryLength() const {
	return age;
}

double Candidate::getLastStep() const {
	return lastStep;
}

double Candidate::getNextStep() const {
	return nextStep;
}

Candidate::Status Candidate::getStatus() const {
	return status;
}

void Candidate::setLastStep(double lstep) {
	lastStep = lstep;
}

void Candidate::setNextStep(double nstep) {
	nextStep = nstep;
}

void Candidate::setStatus(Status stat) {
	status = stat;
}

} // namespace mpc

#endif /* CANDIDATE_H_ */

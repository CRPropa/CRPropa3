#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include "mpc/ParticleState.h"
#include <vector>
#include <map>
#include <sstream>

namespace mpc {

/**
 @class InteractionState
 @brief Candidate state for stochastic interactions.
 */
struct InteractionState {
	double distance;
	int channel;
};

/**
 @class Candidate
 @brief All information about the the cosmic ray candidate.
 */
class Candidate {
public:
	enum Status {
		Active = 0, Detected, OutOfBounds, Stopped, UserDefined
	};

	ParticleState current;
	ParticleState last;
	ParticleState initial;
	std::vector<Candidate *> secondaries;
	std::string history;

private:
	double redshift, trajectoryLength;
	double currentStep, nextStep;
	Status status;
	std::map<std::string, InteractionState> interactionStates;

public:
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

	bool getInteractionState(std::string moduleName, InteractionState &state);
	void setInteractionState(std::string moduleName, InteractionState state);
	const std::map<std::string, InteractionState> getInteractionStates() const;
	void clearInteractionStates();

	void addSecondary(int id, double energy);
};

} // namespace mpc

#endif /* CANDIDATE_H_ */

#ifndef MPC_CANDIDATE_H_
#define MPC_CANDIDATE_H_

#include "mpc/ParticleState.h"
#include "mpc/Referenced.h"

#include <vector>
#include <map>
#include <sstream>

namespace mpc {

/**
 @class InteractionState
 @brief Candidate state for stochastic interactions.
 */
struct InteractionState {
	InteractionState() :
			distance(0), channel(0) {

	}
	InteractionState(double distance, int channel) :
			distance(distance), channel(channel) {

	}
	double distance;
	int channel;
};

/**
 @class Candidate
 @brief All information about the the cosmic ray candidate.
 */
class Candidate: public Referenced {
public:
	bool active;
	ParticleState current;
	ParticleState initial;
	std::vector<ref_ptr<Candidate> > secondaries;

private:
	double redshift, trajectoryLength;
	double currentStep, nextStep;
	std::map<std::string, std::string> properties;
	std::map<std::string, InteractionState> interactionStates;

public:
	Candidate();
	Candidate(const ParticleState &state);
	virtual ~Candidate() {};

	bool isActive() const;
	void setActive(const bool b);

	double getRedshift() const;
	void setRedshift(double z);

	double getTrajectoryLength() const;
	void setTrajectoryLength(double a);

	double getCurrentStep() const;
	void setCurrentStep(double lstep);

	double getNextStep() const;
	void setNextStep(double step);
	void limitNextStep(double step);

	bool getProperty(const std::string &name, std::string &value) const;
	void setProperty(const std::string &name, const std::string &value);
	bool hasProperty(const std::string &name) const;

	bool getInteractionState(const std::string &moduleName,
			InteractionState &state) const;
	void setInteractionState(const std::string &moduleName,
			const InteractionState &state);
	const std::map<std::string, InteractionState> getInteractionStates() const;
	void clearInteractionStates();

	void addSecondary(int id, double energy);
	void clearSecondaries();
};

} // namespace mpc

#endif /* MPC_CANDIDATE_H_ */

#ifndef MPC_CANDIDATE_H_
#define MPC_CANDIDATE_H_

#include "mpc/ParticleState.h"
#include "mpc/Referenced.h"
#include "mpc/AssocVector.h"

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
	double distance; /**< Free distance of the interaction in [m] comoving units */
	int channel; /**< Interaction ID */
};

/**
 @class Candidate
 @brief All information about the cosmic ray.

 The Candidate is a passive object, that holds all information about the initial and current state of the cosmic ray and the status of propagation.
 */
class Candidate: public Referenced {
public:
	ParticleState initial; /**< Particle state at the beginning of propagation */
	ParticleState current; /**< Current particle state */
	ParticleState previous; /**< Particle state at the end of the previous step */

	std::vector<ref_ptr<Candidate> > secondaries; /**< Secondary particles created in interactions */

	typedef Loki::AssocVector<std::string, std::string> PropertyMap;
	typedef Loki::AssocVector<std::string, InteractionState> InteractionStatesMap;

private:
	bool active; /**< Active status */
	double redshift; /**< Current simulation time-point in terms of redshift z, z = 0 being the present */
	double trajectoryLength; /**< Comoving distance [m] the candidate has travelled so far */
	double currentStep; /**< Size of the currently performed step in [m] comoving units */
	double nextStep; /**< Proposed size of the next propagation step in [m] comoving units */

	PropertyMap properties; /**< Map of property names and their values. */
	InteractionStatesMap interactionStates; /**< Map of physical interactions that are scheduled to happen to the candidate. */

public:
	Candidate();
	/** Creates a candidate, initializing the initial, previous and current particle state with the argument. */
	Candidate(const ParticleState &state);

	bool isActive() const;
	void setActive(bool b);

	void setTrajectoryLength(double length);
	double getTrajectoryLength() const;

	void setRedshift(double z);
	double getRedshift() const;

	/** Sets the current step and increases the trajectory length accordingly. Only the propagation module should use this. */
	void setCurrentStep(double step);
	double getCurrentStep() const;

	/** Sets the proposed next step. Only the propagation module should use this. */
	void setNextStep(double step);
	double getNextStep() const;
	void limitNextStep(double step);

	void setProperty(const std::string &name, const std::string &value);
	bool removeProperty(const std::string &name);
	bool getProperty(const std::string &name, std::string &value) const;
	bool hasProperty(const std::string &name) const;

	bool getInteractionState(const std::string &moduleName,
			InteractionState &state) const;
	void setInteractionState(const std::string &moduleName,
			const InteractionState &state);
	const InteractionStatesMap getInteractionStates() const;
	void clearInteractionStates();

	void addSecondary(int id, double energy);
	void clearSecondaries();
};

} // namespace mpc

#endif /* MPC_CANDIDATE_H_ */

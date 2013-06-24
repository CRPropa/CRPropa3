#ifndef CRPROPA_CANDIDATE_H
#define CRPROPA_CANDIDATE_H

#include "crpropa/ParticleState.h"
#include "crpropa/Referenced.h"
#include "crpropa/AssocVector.h"

#include <vector>
#include <map>
#include <sstream>

namespace crpropa {

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

 The Candidate is a passive object, that holds the information about the state
 of the cosmic ray at the beginning of propagation and in the current and
 previous propagation step.
 It holds status information about the cosmic ray and the simulation itself.
 */
class Candidate: public Referenced {
public:
	ParticleState source; /**< Particle state at the source */
	ParticleState created; /**< Particle state at the creation */
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
	InteractionStatesMap interactionStates; /**< Map of interactions that are scheduled to happen to the candidate. */

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

	bool getInteractionState(const std::string &name,
			InteractionState &state) const;
	void setInteractionState(const std::string &name,
			const InteractionState &state);
	const InteractionStatesMap getInteractionStates() const;
	void removeInteractionState(const std::string &name);
	void clearInteractionStates();

	void addSecondary(int id, double energy);
	void clearSecondaries();
};

} // namespace crpropa

#endif // CRPROPA_CANDIDATE_H

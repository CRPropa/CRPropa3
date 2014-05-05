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
 @class Candidate
 @brief All information about the cosmic ray.

 The Candidate is a passive object, that holds the information about the state
 of the cosmic ray and the simulation itself.
 */
class Candidate: public Referenced {
public:
	ParticleState source; /**< Particle state at the source */
	ParticleState created; /**< Particle state of parent particle at the time of creation */
	ParticleState current; /**< Current particle state */
	ParticleState previous; /**< Particle state at the end of the previous step */

	std::vector<ref_ptr<Candidate> > secondaries; /**< Secondary particles from interactions */

	typedef Loki::AssocVector<std::string, std::string> PropertyMap;
	PropertyMap properties; /**< Map of property names and their values. */

private:
	bool active; /**< Active status */
	double redshift; /**< Current simulation time-point in terms of redshift z */
	double trajectoryLength; /**< Comoving distance [m] the candidate has travelled so far */
	double currentStep; /**< Size of the currently performed step in [m] comoving units */
	double nextStep; /**< Proposed size of the next propagation step in [m] comoving units */

public:
	Candidate();
	Candidate(const ParticleState &state); /**< Creates a candidate, initializing the initial, previous and current particle state with the argument. */

	bool isActive() const;
	void setActive(bool b);

	void setTrajectoryLength(double length);
	double getTrajectoryLength() const;

	void setRedshift(double z);
	double getRedshift() const;

	void setCurrentStep(double step); /**< Sets the current step and increases the trajectory length accordingly. Only the propagation module should use this. */
	double getCurrentStep() const;

	void setNextStep(double step); /**< Sets the proposed next step. Only the propagation module should use this. */
	double getNextStep() const;
	void limitNextStep(double step);

	void setProperty(const std::string &name, const std::string &value);
	bool getProperty(const std::string &name, std::string &value) const;
	bool removeProperty(const std::string &name);
	bool hasProperty(const std::string &name) const;

	void addSecondary(int id, double energy);
	void clearSecondaries();
};

} // namespace crpropa

#endif // CRPROPA_CANDIDATE_H

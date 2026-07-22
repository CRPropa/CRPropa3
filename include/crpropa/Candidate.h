#ifndef CRPROPA_CANDIDATE_H
#define CRPROPA_CANDIDATE_H

#include "crpropa/ParticleState.h"
#include "crpropa/Referenced.h"
#include "crpropa/Variant.h"

#include <unordered_map>
#include <vector>
#include <map>
#include <sstream>
#include <stdint.h>

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class Candidate Candidate.h include/crpropa/Candidate.h
 @brief All information about the cosmic ray.

 The Candidate is a passive object, that holds the information about the state
 of the cosmic ray and the simulation itself.
 */
class Candidate {
public:
	ParticleState source; /**< Particle state at the source */
	ParticleState created; /**< Particle state of parent particle at the time of creation */
	ParticleState current; /**< Current particle state */
	ParticleState previous; /**< Particle state at the end of the previous step */

	std::vector<ref_ptr<Candidate> > secondaries; /**< Secondary particles from interactions */

	typedef std::unordered_map<std::string, Variant> PropertyMap;
	PropertyMap properties; /**< Map of property names and their values. */

	/** Parent candidate. 0 if no parent (initial particle). Must not be a ref_ptr to prevent circular referencing. */
	Candidate *parent;

private:
	bool active; /**< Active status */
	double weight; /**< Weight of the candidate */
	double redshift; /**< Current simulation time-point in terms of redshift z */
	double trajectoryLength; /**< Comoving distance [m] the candidate has traveled so far */
	double currentStep; /**< Size of the currently performed step in [m] comoving units */
	double nextStep; /**< Proposed size of the next propagation step in [m] comoving units */
	std::string tagOrigin; /**< Name of interaction/source process which created this candidate*/
	double time; /**< Time [s] that has passed in the laboratory frame of reference */

	static uint64_t nextSerialNumber;
	uint64_t serialNumber;

public:
	/** Constructor
	 * @param id  Particle ID after the 2012 Monte Carlo nuclear code scheme, use nuleusId(A,Z) for nuclei
	 * @param energy  Particle energy
	 * @param position  Start position
	 * @param direction  Start direction
	 * @param z  Redshift
	 * @param weight  Statistical weight (important property for some modules, usally just 1)
	 * @param tagOrigin  Usually either "PRIM" for primary particle or "SEC" for secondary particle, but can be anything
	 */
	Candidate(
		int id = 0,
		double energy = 0,
		Vector3d position = Vector3d(0, 0, 0),
		Vector3d direction = Vector3d(-1, 0, 0),
		double z = 0,
		double weight = 1., 
		std::string tagOrigin = "PRIM"
	);

	/**
	 Creates a candidate, initializing the Candidate::source, Candidate::created,
	 Candidate::previous and Candidate::current state with the argument.
	 @param state  ParticleState for source, created, previous and current. Makes copies.
	 */
	Candidate(const ParticleState &state);

	/** Checks if particle is still active */
	bool isActive() const;
	/** Sets particle active or unactive
	 * When the particle is set unactive it will first finish its current step,
	 * only at the beginning of the next step the particle will be recognized as unactive
	 * and the simulation for that particle will stop.
	 * It might be possible that some modules additionally check if the particle is still active,
	 * if that happens, it might allready be ignored in those modules during the same step,
	 * this depends on the order the modules are added to ModuleList.
	 * @param b  Activate state of particle, false=deactivated
	 */
	void setActive(bool b);

	/** Sets trajectory length
	 * Mostly used by propagators over setCurrentStep, but usefull if particle should start
	 * at a later position.
	 * @param length  Trajectory length in meter
	 */
	void setTrajectoryLength(double length);
	/** Returns current trajectory length */
	double getTrajectoryLength() const;
	
	/** Returns absolute of current velocity
	 * To get the current velocity vector you can use Candidate.current.getVelocity()
	 */
	double getVelocity() const;

	void setRedshift(double z);
	double getRedshift() const;

	/**
	 Sets weight of each candidate.
	 Weights are calculated for each tracked secondary.
	 */
	void setWeight(double weight);
	/** Updates Weight
	 * Multiplies the current weight with the given weight
	 */
    void updateWeight(double weight);
	double getWeight() const;

	/**
	 Sets the current step and increases the trajectory length and time accordingly.
	 Only the propagation module should use this.
	 @param step  Current step in meter
	 */
	void setCurrentStep(double step);
	/** @return Current stepsize in meter */
	double getCurrentStep() const;

	/**
	 Sets the proposed next step.
	 Only the propagation module should use this.
	 @param step  Proposed next stepsize in meter
	 */
	void setNextStep(double step);
	/** @return  Proposed next stepsize in meter */
	double getNextStep() const;

	/**
	 Sets the tagOrigin of the candidate. Can be used to trace back the interactions
	 */
	void setTagOrigin(std::string tagOrigin);
	std::string getTagOrigin() const;

	/**
	 Sets the time of the candidate.
	 This is done automatically together with the increase of TrajectoryLength,
	 since CRPropa assumes lightspeed in every case both TrajectoryLength and Time are equal.
	 @param t  Time in seconds
	 */
	void setTime(double t);
	/** Returns the time of the candidate.
	 * The time is tracked alongside TrajectoryLength by dividing the current TrajecoryLength by c
	 * @return Current time in seconds
	 */
	double getTime() const;

	/**
	 Make a bid for the next step size: the lowest wins.
	 @param step  The bid in meter
	 */
	void limitNextStep(double step);

	/** Sets a arbitrary property
	 * This function either creates a property if it does not exist or updates it.
	 * @param key  Key to put into unordered_map
	 * @param value  Any Variant object, Variant can represent a variety of data types so
	 * that a property can have any basic datatype
	 */
	void setProperty(const std::string &key, const Variant &value);
	/** Returns the value of the Property with given key
	 * This function loops through the unordered_map and returns the first property fitting to
	 * the provided key. If no key value pair is found it throws an error
	 * @param key  Key to search for
	 */
	const Variant &getProperty(const std::string &key) const;
	/** Tries to remove property that has given key
	 * @return Returns true if property was removed successfully or false if not
	 */
	bool removeProperty(const std::string &key);
	/** @return Returns true if property exists, false otherwise */
	bool hasProperty(const std::string &key) const;

	/** Add a new candidate to the list of secondaries.
	 Adds a new candidate to the list of secondaries of this candidate.
	 The secondaries Candidate::source and Candidate::previous state are set to the _source_ and _previous_ state of its parent.
	 The secondaries Candidate::created and Candidate::current state are set to the _current_ state of its parent, except for the secondaries current energy and particle id.
	 Trajectory length and redshift are copied from the parent.
	 @param c Candidate
	 */
	void addSecondary(ref_ptr<Candidate> c);
	/**
	 Add a new candidate to the list of secondaries.
	 @param id			particle ID of the secondary after the 2012 Monte Carlo nuclear code scheme
	 @param energy		energy of the secondary
	 @param w			weight of the secondary
	 @param tagOrigin 	tag of the secondary
	 */
	void addSecondary(int id, double energy, double w = 1., std::string tagOrigin = "SEC");
	/**
	 Add a new candidate to the list of secondaries.
	 @param id			particle ID of the secondary after the 2012 Monte Carlo nuclear code scheme
	 @param energy		energy of the secondary
	 @param position	start position of the secondary
	 @param w			weight of the secondary
	 @param tagOrigin 	tag of the secondary
	 */
	void addSecondary(int id, double energy, Vector3d position, double w = 1., std::string tagOrigin = "SEC");
	/** Clears all stored secondaries (deletes them) */
	void clearSecondaries();

	std::string getDescription() const;

	/** @return Returns a the unique serial number of the particle */
	uint64_t getSerialNumber() const;
	/** Sets a custom serial number
	 * @param snr  Custom serial number
	 */
	void setSerialNumber(const uint64_t snr);

	/** @return Serial number of candidate at source*/
	uint64_t getSourceSerialNumber() const;

	/** @return Serial number of candidate at creation */
	uint64_t getCreatedSerialNumber() const;

	/** Set the next serial number to use */
	static void setNextSerialNumber(uint64_t snr);

	/** @return Get the next serial number that will be assigned */
	static uint64_t getNextSerialNumber();

	/** Create an exact clone of candidate
	 @param recursive	recursively clone and add the secondaries
	 */
	ref_ptr<Candidate> clone(bool recursive = false) const;

	/**
	 Copy the source particle state to the current state
	 and activate it if inactive, e.g. restart it
	*/
	void restart();
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_CANDIDATE_H

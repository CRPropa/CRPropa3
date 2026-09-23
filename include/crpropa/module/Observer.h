#ifndef CRPROPA_OBSERVER_H
#define CRPROPA_OBSERVER_H

#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "../Candidate.h"
#include "../Module.h"
#include "../Referenced.h"
#include "../Vector3.h"
#include "../Geometry.h"

namespace crpropa {

enum DetectionState {
	DETECTED, VETO, NOTHING
};

/** \addtogroup Observer
 * @{
 */


/**
 @class ObserverFeature
 @brief Abstract base class for features of observers
 */
class ObserverFeature: public Referenced {
protected:
	std::string description;
public:
	virtual DetectionState checkDetection(Candidate *candidate) const;
	virtual void onDetection(Candidate *candidate) const;
	virtual std::string getDescription() const;
};


/**
 @class Observer
 @brief General particle observer

 In a UHECR propagation simulation, particles are stepped forward through a
 (potentially enormous) simulated volume, but only a small fraction of that
 trajectory is of physical interest -- for instance, the moment a particle
 crosses a boundary representing an observer's position, or enters a
 specific region of interest. The Observer class implements this selection:
 attached to the simulation, it checks at every propagation step whether
 the current particle satisfies one or more detection conditions, similar
 in spirit to a detector in a physical experiment recording only the
 particles that actually reach it.

 Multiple ObserverFeature conditions can be attached to one Observer, each
 independently classifying a particle as DETECTED, VETO, or having no
 opinion (NOTHING) -- for example, ObserverSurface for a bounding geometry,
 or veto-type features that suppress detection under certain conditions
 (e.g. an inactive candidate or a specific particle type). A single VETO
 from any feature overrides a DETECTED result from any other, regardless
 of the order in which features are checked, so vetoes act as a safeguard
 against reporting particles that satisfy a detection condition but
 shouldn't actually count as observed. Only when the combined result is
 DETECTED are any attached actions triggered: each feature's own
 onDetection() hook runs, an optional user-supplied output module can
 process the candidate (optionally on a cloned copy, leaving the original
 still propagating), a property flag can be set on the candidate, and --
 if configured via setDeactivateOnDetection() -- the candidate can be
 deactivated so it is no longer propagated further, the standard way of
 stopping a particle's simulation once it has effectively "arrived" at the
 observer.

 Because a DETECTED result from a single feature is enough to trigger the
 observer's actions -- provided no other feature returns VETO -- different
 ObserverFeature types are often combined on one Observer to express a
 compound condition (e.g. a detection surface plus a veto that excludes
 certain particle types). For genuinely independent observation goals,
 separate Observer instances are typically used instead -- for example,
 one Observer using ObserverTracking to record a particle's full
 trajectory through a region, alongside a separate Observer using
 ObserverSurface to stop propagation once the particle exits the
 simulated volume.

 @see ObserverFeature, ObserverSurface, ObserverTracking, Source
 */
class Observer: public Module {
	std::string flagKey;
	std::string flagValue;
private:
	std::vector<ref_ptr<ObserverFeature> > features;
	ref_ptr<Module> detectionAction;
	bool clone;
	bool makeInactive;
public:
	/** Default observer constructor
	 */
	Observer();
	/** Add a feature to the observer
	 @param feature		observer feature to be added to the Observer object
	 */
	void add(ObserverFeature *feature);
	/** Perform some specific actions upon detection of candidate
	 @param action		module that performs a given action when candidate is detected
	 @param clone		if true, clone candidate
	 */
	void onDetection(Module *action, bool clone = false);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void setFlag(std::string key, std::string value);
	/** Determine whether candidate should be deactivated on detection
	 @param deactivate	if true, deactivate detected particles; if false, continue tracking them
	 */
	void setDeactivateOnDetection(bool deactivate);
};


/**
 @class ObserverDetectAll
 @brief Detects all particles
 */
class ObserverDetectAll: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverSurface
 @brief Detects particles crossing the boundaries of a defined surface (see, e.g., `Geometry` module)
 */
class ObserverSurface: public ObserverFeature {
private:
	ref_ptr<Surface> surface;
public:
	/** Constructor
	 @param surface		object with some specific geometric (see Geometry.h)
	*/
	ObserverSurface(Surface* surface);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverTracking
 @brief Tracks particles inside a sphere
 */
class ObserverTracking: public ObserverFeature {
private:
	Vector3d center;
	double radius;
    double stepSize;
public:
	/** Constructor
	 @param center		vector containing the coordinates of the center of the sphere
	 @param radius		radius of the sphere
	 @param stepSize	observer will keep track of particles at every step with this size
	*/
	ObserverTracking(Vector3d center, double radius, double stepSize = 0);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class Observer1D
 @brief Detects particles when reaching x = 0

 This module detects particles when reaching x = 0 and also limits the next step size to prevent candidates from overshooting.
 */
class Observer1D: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverRedshiftWindow
 @brief Detects particles in a given redshift window

 When added to an observer, this feature generalizes it to four dimensions.
 The fourth dimension is the redshift, a proxy for time. This is particularly
 useful in "4D" studies, including either time-dependence (e.g. flaring objects),
 or in 3D studies including cosmological evolution.
 Note that redshifts should be assigned to sources when using this feature.
 This can be done with: SourceRedshift, SourceRedshift1D, SourceUniformRedshift,
 and SourceRedshiftEvolution.
 */
class ObserverRedshiftWindow: public ObserverFeature {
private:
	double zmin, zmax;
public:
	/** Constructor
	 @param zmin	lower bound of redshift interval
	 @param zmax	upper bound of redshift interval
	 */
	ObserverRedshiftWindow(double zmin = 0, double zmax = 0.1);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverInactiveVeto
 @brief Veto for inactive candidates
 */
class ObserverInactiveVeto: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverNucleusVeto
 @brief Veto for nuclei (including protons and neutrons)
 */
class ObserverNucleusVeto: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverNeutrinoVeto
 @brief Veto for neutrinos
 */
class ObserverNeutrinoVeto: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverPhotonVeto
 @brief Veto for photons
 */
class ObserverPhotonVeto: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverElectronVeto
 @brief Veto for electrons and positrons
 */
class ObserverElectronVeto: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverParticleIdVeto
 @brief Custom veto for user-defined particle types
 Vetoes for more than one type of particle can be added by calling this
 feature multiple times.
 */
class ObserverParticleIdVeto: public ObserverFeature {
private:
	int vetoParticleId;
public:
	/** Constructor
	 @param id		id of the particle following the PDG numbering scheme
	 */
	ObserverParticleIdVeto(int id);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverTimeEvolution
 @brief Observes the time evolution of the candidates (phase-space elements)
 This observer is very useful if the time evolution of the particle density is needed. It detects all candidates in lin-spaced, log-spaced, or user-defined time intervals and limits the nextStep of candidates to prevent overshooting of detection intervals.
 */
class ObserverTimeEvolution: public ObserverFeature {
protected:
	int nIntervals;  // number of time invervals
	bool isLogarithmicScaling = false;  // enables or disables logarithmic scaling for the intervals
	bool doDetListConstruction = true;  // enables the construction of detList in the relevant functions (addTime, addTimeRange)
	double minimum;  // the minimum time
	double maximum;  // the maximum time
	/** Vector containing all used times. 
	 It is only constructed by the user manually.
	 If it is not empty, the vector will be used instead of the getTime function.
	 (leave empty if you want to rather use functions)
	*/
	std::vector<double> detList;
	/**
	 A temporary storage for detList, this enables the return of a List in getTimes
	 without risking to modify detList
	 */
	mutable std::vector<double> tempDetList;
	
public:
	/** Default constructor
	 */
	ObserverTimeEvolution();
	/** Constructor
	 @param min		minimum time
	 @param dist	time interval for detection
	 @param numb	number of time intervals

	 This constructor calculates the maximum from max = min + (numb - 1) * dist
	 */
	ObserverTimeEvolution(double min, double dist, double numb);
	/** Constructor
	 @param min		minimum time
	 @param max	    maximum time
	 @param numb	number of time intervals
	 @param log     log (input: true) or lin (input: false) scaling between min and max with numb steps
	 
	 This constructor sets the maximum directly and gets numb automatically.
	 You need to set the log parameter, since an overload for the first three doubles exist.
	 */
	ObserverTimeEvolution(double min, double max, double numb, bool log);
	/** Constructor
	 @param detList	user defined vector<double> with times to check

	 This constructor uses a predefined vector containing the times that should be observed.
	 The so created detList can then be modified via addTime, addTimeRange and setTimes.
	 */
	ObserverTimeEvolution(const std::vector<double> &detList);
	/** Destructor
	 */
	~ObserverTimeEvolution(){}

	/** Function
	 Generates the detList if it is empty when for example the 
	 ObserverTimeEvolution(const std::vector<double> &detList) constructor
	 was not used.
	 Use this function to create a detList with can then be modified.
	 When detList is not empty its entries are used instead of a runtime calculation of the times.	 
	 */
	void constructDetListIfEmpty();
	/** Function
	 @param candidate	Candidate usally given by a module list

	 Checks whether to make a detection at the current step of candidate or not.
	 This function is called in Observer.process with the simulated Candidate.
	 */
	DetectionState checkDetection(Candidate *candidate) const;
	/** Function
	 @param enableConstruction	if true, constructs detList from range of min, max, numb
	 when calling addTime
	 Clears the content of detList
	 */
	void clear();
	/** Function
	 Checks if detList is empty
	 */
	bool empty(){return detList.empty();}

	// setter functions:
	/** Function
	 @param time	Time that should be appended to detList

	 Makes a push_back on detList with the given time.	 
	 */
	void addTime(const double &time);
	/** Function
	 @param min		minimum time
	 @param max	    maximum time
	 @param numb	number of time intervals
	 @param log     log (input: true) or lin (input: false) scaling between min and max with numb steps

	 Appends a linear or logarithmic time range to detList via repeatedly calling addTime
	 */
	void addTimeRange(double min, double max, double numb, bool log = false);
	/** Function
	 @param detList	vector<double> containing times when to observe

	 Sets this->detList to detList, with this it is possible to fully modify detList
	 */
	void setTimes(const std::vector<double> &detList);
	void setMinimum(double min);
	void setMaximum(double max){this->maximum = max;}
	void setNIntervals(int numb){this->nIntervals = numb;}
	void setIsLogarithmicScaling(bool log){this->isLogarithmicScaling = log;}

	// getter functions:
	/** Function
	 @param index	index of the required time

	 Replaces the previous preconstructed detList and return the time for a specific index
	 */
	virtual double getTime(std::size_t index) const;
	double getMinimum() const {return minimum;}
	double getMaximum() const {return maximum;}
	int getNIntervals() const {return nIntervals;}
	bool getIsLogarithmicScaling() const {return isLogarithmicScaling;}
	/** Function 
	 Returns a vector<double> containing all times between min and max generated by getTime.
	 This function does not return detList directly, it rather appends a new vector with
	 getTime numb times.
	 */
	const std::vector<double>& getTimes() const;
	/** Function
	 Returns a string containing a representation of all times.
	 This function does not create a detList.
	 */
	std::string getDescription() const;
};

/** @} */

}

#endif // CRPROPA_OBSERVER_H

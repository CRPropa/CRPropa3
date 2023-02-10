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
 @class ObserverPoint
 @brief Detects particles when reaching x = 0

Should be removed and replaced by Observer1D
 */
class ObserverPoint: public ObserverFeature {
public:
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
private:
	std::vector<double> detList;
public:
	/** Default constructor
	 */
	ObserverTimeEvolution();
	/** Constructor
	 @param min		minimum time
	 @param dist	time interval for detection
	 @param numb	number of time intervals
	 */
	ObserverTimeEvolution(double min, double dist, double numb);
	/** Constructor
	 @param min		minimum time
	 @param max	    maximum time
	 @param numb	number of time intervals
	 @param log     log (input: true) or lin (input: false) scaling between min and max with numb steps
	 */
	ObserverTimeEvolution(double min, double max, double numb, bool log);
	// Add a new time step to the detection time list of the observer
	void addTime(const double &position);
	// Using log or lin spacing of times in the range between min and
	// max for observing particles
	void addTimeRange(double min, double max, double numb, bool log = false);
	const std::vector<double>& getTimes() const;
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};
/** @} */

}

#endif // CRPROPA_OBSERVER_H

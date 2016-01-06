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

namespace crpropa {

enum DetectionState {
	DETECTED, VETO, NOTHING
};

/**
 @class ObserverFeature
 @brief Abstract base class for features of cosmic ray observers
 */
class ObserverFeature: public Referenced {
protected:
	std::string description;
public:
	virtual DetectionState checkDetection(Candidate *candidate) const;
	virtual void onDetection(Candidate *candidate) const;
	virtual std::string getDescription() const;
	virtual void beginRun();
	virtual void endRun();
};

/**
 @class Observer
 @brief General cosmic ray observer
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
	Observer();
	void add(ObserverFeature *feature);
	void onDetection(Module *action, bool clone = false);
	void beginRun();
	void endRun();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void setFlag(std::string key, std::string value);
	void setDeactivateOnDetection(bool deactivate);
};

/**
 @class ObserverSmallSphere
 @brief Detects particles upon entering a sphere
 */
class ObserverSmallSphere: public ObserverFeature {
private:
	Vector3d center;
	double radius;
public:
	ObserverSmallSphere(Vector3d center = Vector3d(0.), double radius = 0);
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
	ObserverTracking(Vector3d center, double radius, double stepSize = 0);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ObserverLargeSphere
 @brief Detects particles upon exiting a sphere
 */
class ObserverLargeSphere: public ObserverFeature {
private:
	Vector3d center;
	double radius;
public:
	ObserverLargeSphere(Vector3d center = Vector3d(0.), double radius = 0);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ObserverPoint
 @brief Detects particles when reaching x = 0

 This module limits the next step size to prevent candidates from overshooting.
 Should be renamed to Observer1D, once old observer-scheme is removed.
 */
class ObserverPoint: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ObserverRedshiftWindow
 @brief Detects particles in a given redshift window
 */
class ObserverRedshiftWindow: public ObserverFeature {
private:
	double zmin, zmax;
public:
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


////////////////////////////////////////////////////////////////////////////////
// old observer scheme

class SmallObserverSphere: public Module {
private:
	Vector3d center;
	double radius;
	std::string flag;
	std::string flagValue;
	bool makeInactive;
	double maximumTrajectory;

public:
	SmallObserverSphere(Vector3d center = Vector3d(0.), double radius = 0,
			std::string flag = "Detected", std::string flagValue = "",
			bool makeInactive = true);
	void process(Candidate *candidate) const;
	void setCenter(Vector3d center);
	void setRadius(double radius);
	void setFlag(std::string flag, std::string flagValue);
	void setMakeInactive(bool makeInactive);
	void setMaximumTrajectory(double maximumTrajectory);
	std::string getDescription() const;
};

/**
 @class LargeObserverSphere
 @brief Flags particles when leaving the sphere.

 Particles are detected when the current position is outside and the previous position inside the sphere.
 In this case the candidate is by default flagged "Detected" made inactive.
 This module limits the next step size to prevent candidates from overshooting.
 */
class LargeObserverSphere: public Module {
private:
	Vector3d center;
	double radius;
	std::string flag;
	std::string flagValue;
	bool makeInactive;

public:
	LargeObserverSphere(Vector3d center = Vector3d(0.), double radius = 0,
			std::string flag = "Detected", std::string flagValue = "",
			bool makeInactive = true);
	void process(Candidate *candidate) const;
	void setCenter(Vector3d center);
	void setRadius(double radius);
	void setFlag(std::string flag, std::string flagValue);
	void setMakeInactive(bool makeInactive);
	std::string getDescription() const;
};

/**
 @class OneDimObserver
 @brief Observer for 1D simulations

 Particles are detected once their x-position gets smaller than a 0.
 In this case the candidate is by flagged "Detected" made inactive.
 This module limits the next step size to prevent candidates from overshooting.
 */
class Observer1D: public Module {
public:
	Observer1D();
	void process(Candidate *candidate) const;
};

/**
 @class DetectAll
 @brief Flags all

 Module that flags every particle it processes.
 Optionally inactivates the particle.
 */
class DetectAll: public Module {
private:
	std::string flag;
	std::string flagValue;
	bool makeInactive;

public:
	DetectAll(std::string flag = "Detected", std::string flagValue = "",
			bool makeInactive = true);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_OBSERVER_H

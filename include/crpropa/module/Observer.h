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
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void setFlag(std::string key, std::string value);
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
	void setCenter(const Vector3d &center);
	void setRadius(float radius);
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

}

#endif // CRPROPA_OBSERVER_H

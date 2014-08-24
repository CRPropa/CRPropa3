#ifndef CRPROPA_OBSERVER_H
#define CRPROPA_OBSERVER_H

#include "crpropa/Module.h"

#include <vector>
#include <fstream>
#include <limits>

namespace crpropa {

/**
 @class ObserverFeature
 @brief Abstract base class for features of cosmic ray observers
 */
class ObserverFeature: public Referenced {
public:
	virtual bool process(Candidate &candidate) const = 0;
};

/**
 @class Observer
 @brief General cosmic ray observer
 */
class Observer: public Module {
private:
	std::vector<ref_ptr<ObserverFeature> > features;
	bool makeInactive;
public:
	Observer(bool makeInactive = false) :
			makeInactive(makeInactive) {
	}
	void add(ObserverFeature *property);
	void process(Candidate *candidate) const;
};

/**
 @class ObserverSmallSphere
 @brief Detects particles upon entering a sphere
 */
class ObserverSmallSphere: public ObserverFeature {
private:
	Vector3d center;
	double radius;
	double maximumTrajectory;
public:
	ObserverSmallSphere(Vector3d center = Vector3d(0.), double radius = 0,
			double maximumTrajectory = std::numeric_limits<double>::max());
	bool process(Candidate &candidate) const;
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
	bool process(Candidate &candidate) const;
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
	bool process(Candidate &candidate) const;
};

/**
 @class ObserverPoint
 @brief Detects particles when reaching x = 0

 Should be renamed to Observer1D
 */
class ObserverPoint: public ObserverFeature {
public:
	bool process(Candidate &candidate) const;
};

/**
 @class ObserverOutput3D
 @brief Plain text output of 3D properties
 */
class ObserverOutput3D: public ObserverFeature {
private:
	mutable std::ofstream fout;
	bool legacy;  // toggle CRPropa 2 output format
public:
	ObserverOutput3D(std::string filename, bool legacy = false);
	~ObserverOutput3D();
	bool process(Candidate &candidate) const;
};

/**
 @class ObserverOutput1D
 @brief Plain text output of 1D properties
 */
class ObserverOutput1D: public ObserverFeature {
private:
	mutable std::ofstream fout;
	bool legacy;  // toggle CRPropa 2 output format
public:
	ObserverOutput1D(std::string filename, bool legacy = false);
	~ObserverOutput1D();
	bool process(Candidate &candidate) const;
};

////////////////////////////////////////////////////////////////////////////////

class SmallObserverSphere: public Module {
private:
	Vector3d center;
	double radius;
	std::string flag;
	std::string flagValue;
	bool makeInactive;

public:
	SmallObserverSphere(Vector3d center = Vector3d(0.), double radius = 0,
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
	void updateDescription();

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
 @class Observer1D
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

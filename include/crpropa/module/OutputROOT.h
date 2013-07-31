#ifndef CRPROPA_OUTPUTROOT_H
#define CRPROPA_OUTPUTROOT_H

#include "crpropa/Module.h"

#ifdef CRPROPA_HAVE_ROOT
#include <TFile.h>
#include <TNtuple.h>
#include <TThread.h>

namespace crpropa {

/**
 @class CRPropa2ROOTEventOutput1D
 @brief Records particles that are inactive and have the property 'Detected' to a ROOT file.
 */
class CRPropa2ROOTEventOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	CRPropa2ROOTEventOutput1D(std::string filename);
	~CRPropa2ROOTEventOutput1D();
	void process(Candidate *candidate) const;
};

/**
 @class ROOTEventOutput1D
 @brief Records particles that are inactive and have the property 'Detected' to a ROOT file.
 */
class ROOTEventOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	ROOTEventOutput1D(std::string filename);
	~ROOTEventOutput1D();
	void process(Candidate *candidate) const;
};


/**
 @class CRPropa2ROOTTrajectoryOutput1D
 @brief Saves trajectories to root file.
 */
class CRPropa2ROOTTrajectoryOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	CRPropa2ROOTTrajectoryOutput1D(std::string filename);
	~CRPropa2ROOTTrajectoryOutput1D();
	void process(Candidate *candidate) const;
};


/**
 @class ROOTTrajectoryOutput1D
 @brief Saves trajectories to root file.
 */
class ROOTTrajectoryOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	ROOTTrajectoryOutput1D(std::string filename);
	~ROOTTrajectoryOutput1D();
	void process(Candidate *candidate) const;
};



/**
 @class CRPropa2ROOTEventOutput3D
 @brief Records particles that have the property 'Detected' to a ROOT file in 3D.
 */
class CRPropa2ROOTEventOutput3D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	CRPropa2ROOTEventOutput3D(std::string filename);
	~CRPropa2ROOTEventOutput3D();
	void process(Candidate *candidate) const;
};

/**
 @class ROOTEventOutput3D
 @brief Records particles that have the property 'Detected' to a ROOT file in 3D.
 */
class ROOTEventOutput3D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	ROOTEventOutput3D(std::string filename);
	~ROOTEventOutput3D();
	void process(Candidate *candidate) const;
};


/**
 @class CRPropa2ROOTTrajectoryOutput3D
 @brief Saves trajectories to root file in 3D.
 */
class CRPropa2ROOTTrajectoryOutput3D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	CRPropa2ROOTTrajectoryOutput3D(std::string filename);
	~CRPropa2ROOTTrajectoryOutput3D();
	void process(Candidate *candidate) const;
};

/**
 @class ROOTTrajectoryOutput3D
 @brief Saves trajectories to root file in 3D.
 */
class ROOTTrajectoryOutput3D: public Module {
	mutable TFile *ROOTFile;
	mutable TNtuple *Ntuple;

public:
	ROOTTrajectoryOutput3D(std::string filename);
	~ROOTTrajectoryOutput3D();
	void process(Candidate *candidate) const;
};


} // namespace crpropa

#endif // CRPROPA_HAVE_ROOT
#endif // CRPROPA_OUTPUTROOT_H

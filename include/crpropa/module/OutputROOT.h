#ifndef CRPROPA_OUTPUTROOT_H
#define CRPROPA_OUTPUTROOT_H

#include "crpropa/Module.h"

#ifdef CRPROPA_HAVE_ROOT
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TThread.h>

namespace crpropa {


/**
 @class ROOTEventOutput1D
 @brief Records particles that are inactive and have the property 'Detected' to a ROOT file.
 */
class ROOTEventOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TTree *Tree;

	mutable int Particle_Type, Initial_Type;
	mutable float Energy_EeV, Initial_Energy_EeV;
	mutable float TrajectoryLength_Mpc;
public:
	ROOTEventOutput1D(std::string filename);
	~ROOTEventOutput1D();
	void process(Candidate *candidate) const;
	void endRun();
};

/**
 @class ROOTTrajectoryOutput1D
 @brief Saves trajectories to root file.
 */
class ROOTTrajectoryOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TTree *Tree;

	mutable int Particle_Type;
	mutable float Energy_EeV;
	mutable float Position_Mpc;
public:
	ROOTTrajectoryOutput1D(std::string filename);
	~ROOTTrajectoryOutput1D();
	void process(Candidate *candidate) const;
	void endRun();
};

/**
 @class ROOTEventOutput3D
 @brief Records particles that have the property 'Detected' to a ROOT file in 3D.
 */
class ROOTEventOutput3D: public Module {
	mutable TFile *ROOTFile;
	mutable TTree *Tree;

	mutable int Particle_Type, Initial_Type;
	mutable float Momentum_E_EeV, Initial_Momentum_E_EeV;
	mutable float TrajectoryLength_Mpc;

	mutable float Position_X_Mpc, Position_Y_Mpc, Position_Z_Mpc;
	mutable float Initial_Position_X_Mpc, Initial_Position_Y_Mpc, Initial_Position_Z_Mpc;
	mutable float Direction_X_Mpc, Direction_Y_Mpc, Direction_Z_Mpc;

public:
	ROOTEventOutput3D(std::string filename);
	~ROOTEventOutput3D();
	void process(Candidate *candidate) const;
	void endRun();
};

/**
 @class ROOTTrajectoryOutput3D
 @brief Saves trajectories to root file in 3D.
 */
class ROOTTrajectoryOutput3D: public Module {
	mutable TFile *ROOTFile;
	mutable TTree *Tree;

	mutable int Particle_Type;
	mutable float Energy_EeV;
	mutable float TrajectoryLength_Mpc;

	mutable float Position_X_Mpc, Position_Y_Mpc, Position_Z_Mpc;
	mutable float Direction_X_Mpc, Direction_Y_Mpc, Direction_Z_Mpc;

public:
	ROOTTrajectoryOutput3D(std::string filename);
	~ROOTTrajectoryOutput3D();
	void process(Candidate *candidate) const;
	void endRun();
};

/**
 @class ROOTPhotonOutput1D
 @brief Records EM-particles that are inactive and have the property 'Detected' to a ROOT file.
 */
class ROOTPhotonOutput1D: public Module {
	mutable TFile *ROOTFile;
	mutable TTree *Tree;

	mutable int Particle_Type, Initial_Type, Parent_Type;
	mutable float Energy_EeV, Initial_Energy_EeV, Parent_Energy_EeV;
	mutable float ComovingDistance_Mpc;
public:
	ROOTPhotonOutput1D(std::string filename);
	~ROOTPhotonOutput1D();
	void process(Candidate *candidate) const;
	void endRun();
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



} // namespace crpropa

#endif // CRPROPA_HAVE_ROOT
#endif // CRPROPA_OUTPUTROOT_H

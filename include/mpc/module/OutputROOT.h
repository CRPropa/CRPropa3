#ifndef MPC_OUTPUTROOT_H_
#define MPC_OUTPUTROOT_H_

#include "mpc/Module.h"

#ifdef MPC_HAVE_ROOT
#include <TFile.h>
#include <TNtuple.h>

namespace mpc {

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

} // namespace mpc

#endif // MPC_HAVE_ROOT
#endif // MPC_OUTPUTROOT_H_

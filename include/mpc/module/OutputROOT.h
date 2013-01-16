#ifndef OUTPUTROOT_H_
#define OUTPUTROOT_H_

#include "mpc/Module.h"
#include <fstream>

#ifdef MPC_HAVE_ROOT
#include <TFile.h>
#include <TNtuple.h>

namespace mpc {

/**
 @class EventOutput 1D ROOT
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
 @class TrajectoryOutput 1D ROOT
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
 @class EventOutput 3D ROOT
 @brief Records particles that are inactive and have the property 'Detected' to a ROOT file in 3D.
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
 @class TrajectoryOutput 3D ROOT
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
#endif /* OUTPUTROOT_H_ */

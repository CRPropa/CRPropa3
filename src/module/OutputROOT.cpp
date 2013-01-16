#include "mpc/module/OutputROOT.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>


#ifdef MPC_HAVE_ROOT

namespace mpc {

/////////////////////// ROOT EVENT OUTPUT 1D //////////////////////////////////
ROOTEventOutput1D::ROOTEventOutput1D(std::string filename) {
  ROOTFile = new TFile(filename.c_str(),"RECREATE","CRPropa output data file");
  Ntuple = new TNtuple("events","CRPropa 1D events","Particle_Type:Initial_Type:Initial_Position_Mpc:Inital_Redshift:Initial_Energy_EeV:Time_Mpc:Energy_EeV");
}

ROOTEventOutput1D::~ROOTEventOutput1D() {
  ROOTFile->Write() ;
  ROOTFile->Close() ;
}

void ROOTEventOutput1D::process(Candidate *c) const {
	if (c->isActive())
		return;
	if (c->hasProperty("Detected")) {
	  #pragma omp critical
	  {
	    Ntuple->Fill(c->current.getId(), 
	                 c->initial.getId(),
			 c->initial.getPosition().getX() / Mpc,
			 c->getRedshift(), 
			 c->initial.getEnergy() / EeV,
			 c->getTrajectoryLength() / Mpc,
			 c->current.getEnergy() / EeV);
	  }
	}
}
////////////////////////////////////////////////////////////////////////////////



/////////////////////// ROOT TRAJECTORY OUTPUT 1D //////////////////////////////
ROOTTrajectoryOutput1D::ROOTTrajectoryOutput1D(std::string filename) {
  ROOTFile = new TFile(filename.c_str(),"RECREATE","CRPropa output data file");
  //  Ntuple = new TNtuple("traj","CRPropa 1D trajectories","Particle_Type:Initial_Type:Time_Mpc:Position_Mpc:Energy_EeV");
  Ntuple = new TNtuple("traj","CRPropa 1D trajectories","Particle_Type:Initial_Type:Time_Mpc:Position_Mpc:Energy_EeV");
}

ROOTTrajectoryOutput1D::~ROOTTrajectoryOutput1D() {
  ROOTFile->Write() ;
  ROOTFile->Close() ;
}

void ROOTTrajectoryOutput1D::process(Candidate *c) const {
	  #pragma omp critical
	  {
	    Ntuple->Fill(c->current.getId(), 
	                 c->initial.getId(),
			 c->getTrajectoryLength() / Mpc,
			 c->current.getPosition().getX() / Mpc,
			 c->current.getEnergy() / EeV);
	  }
}  
////////////////////////////////////////////////////////////////////////////////



/////////////////////// ROOT EVENT OUTPUT 3D //////////////////////////////////
ROOTEventOutput3D::ROOTEventOutput3D(std::string filename) {
  ROOTFile = new TFile(filename.c_str(),"RECREATE","CRPropa output data file");
  Ntuple = new TNtuple("events","CRPropa 3D events","Particle_Type:Initial_Type:Initial_Position_X_Mpc:Initial_Position_Y_Mpc:Initial_Position_Z_Mpc:Initial_Momentum_E_EeV:Initial_Momentum_theta:Initial_Momentum_phi:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Momentum_E_EeV:Momentum_theta:Momentum_phi");
}

ROOTEventOutput3D::~ROOTEventOutput3D() {
  ROOTFile->Write() ;
  ROOTFile->Close() ;
}

void ROOTEventOutput3D::process(Candidate *c) const {
	if (c->isActive())
		return;
	if (c->hasProperty("Detected")) {
	  #pragma omp critical
	  {
	    Ntuple->Fill(c->current.getId(), 
	                 c->initial.getId(),
			 c->initial.getPosition().getX() / Mpc,
			 c->initial.getPosition().getY() / Mpc,
			 c->initial.getPosition().getZ() / Mpc,
			 c->initial.getEnergy() / EeV,
			 c->initial.getDirection().getTheta(),
			 c->initial.getDirection().getPhi(),
			 c->getTrajectoryLength() / Mpc,
			 c->current.getPosition().getX() / Mpc,
			 c->current.getPosition().getY() / Mpc,
			 c->current.getPosition().getZ() / Mpc,
			 c->current.getEnergy() / EeV,
			 c->current.getDirection().getTheta(),
			 c->current.getDirection().getPhi());
	  }
	}
}
////////////////////////////////////////////////////////////////////////////////


/////////////////////// ROOT TRAJECTORY OUTPUT 3D //////////////////////////////
ROOTTrajectoryOutput3D::ROOTTrajectoryOutput3D(std::string filename) {
  ROOTFile = new TFile(filename.c_str(),"RECREATE","CRPropa output data file");
  Ntuple = new TNtuple("traj","CRPropa 3D trajectories","Particle_Type:Initial_Type:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Direction_X:Direction_Y:Direction_Z:Energy_EeV");
}

ROOTTrajectoryOutput3D::~ROOTTrajectoryOutput3D() {
  ROOTFile->Write() ;
  ROOTFile->Close() ;
}

void ROOTTrajectoryOutput3D::process(Candidate *c) const {
	  #pragma omp critical
	  {
	    Ntuple->Fill(c->current.getId(), 
	                 c->initial.getId(),
			 c->getTrajectoryLength() / Mpc,
			 c->current.getPosition().getX() / Mpc,
			 c->current.getPosition().getY() / Mpc,
			 c->current.getPosition().getZ() / Mpc,
			 c->current.getDirection().getX(),
			 c->current.getDirection().getY(),
			 c->current.getDirection().getZ(),
			 c->current.getEnergy() / EeV);
	  }
}  
////////////////////////////////////////////////////////////////////////////////

} // namespace mpc
#endif // MPC_HAVE_ROOT

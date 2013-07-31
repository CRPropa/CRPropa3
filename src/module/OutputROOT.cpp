#include "crpropa/module/OutputROOT.h"
#include "crpropa/Cosmology.h"

#ifdef CRPROPA_HAVE_ROOT

namespace crpropa {

/////////////////////// CRPropa2ROOT EVENT OUTPUT 1D //////////////////////////////////
CRPropa2ROOTEventOutput1D::CRPropa2ROOTEventOutput1D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa 2 output data file");
	Ntuple = new TNtuple("events", "CRPropa 1D events",
			     "Particle_Type:Initial_Type:Initial_Position_Mpc:Initial_Energy_EeV:Energy_EeV");
	TThread::UnLock();
}

CRPropa2ROOTEventOutput1D::~CRPropa2ROOTEventOutput1D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void CRPropa2ROOTEventOutput1D::process(Candidate *c) const {
	if (not (c->hasProperty("Detected")))
		return;

	c->removeProperty("Detected");

	TThread::Lock();
#pragma omp critical
	{
	  Ntuple->Fill(convertToCRPropaId(c->current.getId()),
		       convertToCRPropaId(c->source.getId()),
		       comoving2LightTravelDistance(c->source.getPosition().getX() / Mpc),
		       c->source.getEnergy() / EeV,
		       c->current.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////


/////////////////////// CRPropa2ROOT TRAJECTORY OUTPUT 1D //////////////////////////////
CRPropa2ROOTTrajectoryOutput1D::CRPropa2ROOTTrajectoryOutput1D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa 2 output data file");
	Ntuple = new TNtuple("traj", "CRPropa 1D trajectories",
			"Particle_Type:Initial_Type:Position_Mpc:Energy_EeV");
	TThread::UnLock();
}

CRPropa2ROOTTrajectoryOutput1D::~CRPropa2ROOTTrajectoryOutput1D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void CRPropa2ROOTTrajectoryOutput1D::process(Candidate *c) const {
	TThread::Lock();
#pragma omp critical
	{
	  Ntuple->Fill(convertToCRPropaId(c->current.getId()), 
		       convertToCRPropaId(c->source.getId()),
		       comoving2LightTravelDistance(c->current.getPosition().getX() / Mpc),
		       c->current.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

/////////////////////// CRPropa2ROOT EVENT OUTPUT 3D //////////////////////////////////
CRPropa2ROOTEventOutput3D::CRPropa2ROOTEventOutput3D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa 2 output data file");
	Ntuple =
			new TNtuple("events", "CRPropa 3D events",
					"Particle_Type:Initial_Type:Initial_Position_X_Mpc:Initial_Position_Y_Mpc:Initial_Position_Z_Mpc:Initial_Momentum_E_EeV:Initial_Momentum_theta:Initial_Momentum_phi:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Momentum_E_EeV:Momentum_theta:Momentum_phi");
	TThread::UnLock();
}

CRPropa2ROOTEventOutput3D::~CRPropa2ROOTEventOutput3D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void CRPropa2ROOTEventOutput3D::process(Candidate *c) const {
	if (not (c->hasProperty("Detected")))
		return;

	c->removeProperty("Detected");

	TThread::Lock();
#pragma omp critical
	{
	  Ntuple->Fill(convertToCRPropaId(c->current.getId()), 
		       convertToCRPropaId(c->source.getId()),
		       comoving2LightTravelDistance(c->source.getPosition().getX() / Mpc),
		       comoving2LightTravelDistance(c->source.getPosition().getY() / Mpc),
		       comoving2LightTravelDistance(c->source.getPosition().getZ() / Mpc),
		       c->source.getEnergy() / EeV,
		       c->source.getDirection().getTheta(),
		       c->source.getDirection().getPhi(),
		       comoving2LightTravelDistance(c->getTrajectoryLength() / Mpc),
		       comoving2LightTravelDistance(c->current.getPosition().getX() / Mpc),
		       comoving2LightTravelDistance(c->current.getPosition().getY() / Mpc),
		       comoving2LightTravelDistance(c->current.getPosition().getZ() / Mpc),
		       c->current.getEnergy() / EeV,
		       c->current.getDirection().getTheta(),
		       c->current.getDirection().getPhi());
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

/////////////////////// CRPropa2ROOT TRAJECTORY OUTPUT 3D //////////////////////////////
CRPropa2ROOTTrajectoryOutput3D::CRPropa2ROOTTrajectoryOutput3D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa 2 output data file");
	Ntuple =
			new TNtuple("traj", "CRPropa 3D trajectories",
					"Particle_Type:Initial_Type:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Direction_X:Direction_Y:Direction_Z:Energy_EeV");
	TThread::UnLock();
}

CRPropa2ROOTTrajectoryOutput3D::~CRPropa2ROOTTrajectoryOutput3D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void CRPropa2ROOTTrajectoryOutput3D::process(Candidate *c) const {
	TThread::Lock();
#pragma omp critical
	{
		Ntuple->Fill(c->current.getId(),
			     c->source.getId(),
			     comoving2LightTravelDistance(c->getTrajectoryLength() / Mpc),
			     comoving2LightTravelDistance(c->current.getPosition().getX() / Mpc),
			     comoving2LightTravelDistance(c->current.getPosition().getY() / Mpc),
			     comoving2LightTravelDistance(c->current.getPosition().getZ() / Mpc),
			     c->current.getDirection().getX(),
			     c->current.getDirection().getY(),
			     c->current.getDirection().getZ(),
			     c->current.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////





/////////////////////// ROOT EVENT OUTPUT 1D //////////////////////////////////
ROOTEventOutput1D::ROOTEventOutput1D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	Ntuple = new TNtuple("events", "CRPropa 1D events",
	 		     "Particle_Type:Energy_EeV:TrajectoryLength_Mpc:Initial_Type:Initial_Energy_EeV");
	TThread::UnLock();
}


ROOTEventOutput1D::~ROOTEventOutput1D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void ROOTEventOutput1D::process(Candidate *c) const {
	if (not (c->hasProperty("Detected")))
		return;

	c->removeProperty("Detected");

	TThread::Lock();
#pragma omp critical
	{
	Ntuple->Fill(c->current.getId(),
		     c->current.getEnergy() / EeV,
		     c->getTrajectoryLength() / Mpc,
		     c->source.getId(),
		     c->source.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////


/////////////////////// ROOT TRAJECTORY OUTPUT 1D //////////////////////////////
ROOTTrajectoryOutput1D::ROOTTrajectoryOutput1D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	Ntuple = new TNtuple("traj", "CRPropa 1D trajectories",
			     "Position_Mpc:Particle_Type:Energy_EeV");
	TThread::UnLock();
}

ROOTTrajectoryOutput1D::~ROOTTrajectoryOutput1D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void ROOTTrajectoryOutput1D::process(Candidate *c) const {
	TThread::Lock();
#pragma omp critical
	{
	  Ntuple->Fill(c->current.getPosition().getX() / Mpc, 
		       c->current.getId(),
		       c->current.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

/////////////////////// ROOT EVENT OUTPUT 3D //////////////////////////////////
ROOTEventOutput3D::ROOTEventOutput3D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	Ntuple = new TNtuple("events", "CRPropa 3D events",
			     "TrajectoryLength_Mpc:Particle_Type:Initial_Type:Momentum_E_EeV:Initial_Momentum_E_EeV:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Initial_Position_X_Mpc:Initial_Position_Y_Mpc:Initial_Position_Z_Mpc:Direction_X_Mpc:Direction_Y_Mpc:Direction_Z_Mpc:Initial_Direction_X_Mpc:Initial_Direction_Y_Mpc:Initial_Direction_Z_Mpc");
	TThread::UnLock();
}

ROOTEventOutput3D::~ROOTEventOutput3D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void ROOTEventOutput3D::process(Candidate *c) const {
	if (not (c->hasProperty("Detected")))
		return;

	c->removeProperty("Detected");

	TThread::Lock();
#pragma omp critical
	{
	  // Ntuple->Fill(c->getTrajectoryLength() / Mpc,
	  // 	       c->current.getId(),
	  // 	       c->source.getId(),
	  // 	       c->current.getEnergy() / EeV,
	  // 	       c->source.getEnergy() / EeV,
	  // 	       c->current.getPosition().getX() / Mpc,
	  // 	       c->current.getPosition().getY() / Mpc,
	  // 	       c->current.getPosition().getZ() / Mpc,
	  // 	       c->source.getPosition().getX() / Mpc,
	  // 	       c->source.getPosition().getY() / Mpc,
	  // 	       c->source.getPosition().getZ() / Mpc,
	  // 	       c->current.getDirection().getX() / Mpc,
	  // 	       c->current.getDirection().getY() / Mpc,
	  // 	       c->current.getDirection().getZ() / Mpc,
	  // 	       c->source.getDirection().getX() / Mpc,
	  // 	       c->source.getDirection().getY() / Mpc,
	  // 	       c->source.getDirection().getZ() / Mpc);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

/////////////////////// ROOT TRAJECTORY OUTPUT 3D //////////////////////////////
ROOTTrajectoryOutput3D::ROOTTrajectoryOutput3D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	Ntuple = new TNtuple("traj", "CRPropa 3D trajectories",
			     "TrajectoryLength_Mpc:Particle_Type:Energy_EeV:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Direction_X_Mpc:Direction_Y_Mpc:Direction_Z_Mpc");
	TThread::UnLock();
}

ROOTTrajectoryOutput3D::~ROOTTrajectoryOutput3D() {
	TThread::Lock();
	ROOTFile->Write();
	ROOTFile->Close();
	TThread::UnLock();
}

void ROOTTrajectoryOutput3D::process(Candidate *c) const {
	TThread::Lock();
#pragma omp critical
	{
	  Ntuple->Fill(c->getTrajectoryLength() / Mpc,
		       c->current.getId(),
		       c->current.getEnergy() / EeV,
		       c->current.getPosition().getX() / Mpc,
		       c->current.getPosition().getY() / Mpc,
		       c->current.getPosition().getZ() / Mpc,
		       c->current.getDirection().getX(),
		       c->current.getDirection().getY(),
		       c->current.getDirection().getZ());
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////


}// namespace crpropa

#endif // CRPROPA_HAVE_ROOT

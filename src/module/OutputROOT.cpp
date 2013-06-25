#include "crpropa/module/OutputROOT.h"

#ifdef CRPROPA_HAVE_ROOT

namespace crpropa {

/////////////////////// ROOT EVENT OUTPUT 1D //////////////////////////////////
ROOTEventOutput1D::ROOTEventOutput1D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	Ntuple =
			new TNtuple("events", "CRPropa 1D events",
					"Particle_Type:Initial_Type:Initial_Position_Mpc:Inital_Redshift:Initial_Energy_EeV:Time_Mpc:Energy_EeV");
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
		Ntuple->Fill(c->current.getId(), c->source.getId(),
				c->source.getPosition().getX() / Mpc, c->getRedshift(),
				c->source.getEnergy() / EeV, c->getTrajectoryLength() / Mpc,
				c->current.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

/////////////////////// ROOT TRAJECTORY OUTPUT 1D //////////////////////////////
ROOTTrajectoryOutput1D::ROOTTrajectoryOutput1D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	//  Ntuple = new TNtuple("traj","CRPropa 1D trajectories","Particle_Type:Initial_Type:Time_Mpc:Position_Mpc:Energy_EeV");
	Ntuple = new TNtuple("traj", "CRPropa 1D trajectories",
			"Particle_Type:Initial_Type:Time_Mpc:Position_Mpc:Energy_EeV");
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
		Ntuple->Fill(c->current.getId(), c->source.getId(),
				c->getTrajectoryLength() / Mpc,
				c->current.getPosition().getX() / Mpc,
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
	Ntuple =
			new TNtuple("events", "CRPropa 3D events",
					"Particle_Type:Initial_Type:Initial_Position_X_Mpc:Initial_Position_Y_Mpc:Initial_Position_Z_Mpc:Initial_Momentum_E_EeV:Initial_Momentum_theta:Initial_Momentum_phi:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Momentum_E_EeV:Momentum_theta:Momentum_phi");
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
		Ntuple->Fill(c->current.getId(), c->source.getId(),
				c->source.getPosition().getX() / Mpc,
				c->source.getPosition().getY() / Mpc,
				c->source.getPosition().getZ() / Mpc,
				c->source.getEnergy() / EeV,
				c->source.getDirection().getTheta(),
				c->source.getDirection().getPhi(),
				c->getTrajectoryLength() / Mpc,
				c->current.getPosition().getX() / Mpc,
				c->current.getPosition().getY() / Mpc,
				c->current.getPosition().getZ() / Mpc,
				c->current.getEnergy() / EeV,
				c->current.getDirection().getTheta(),
				c->current.getDirection().getPhi());
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

/////////////////////// ROOT TRAJECTORY OUTPUT 3D //////////////////////////////
ROOTTrajectoryOutput3D::ROOTTrajectoryOutput3D(std::string filename) {
	TThread::Lock();
	ROOTFile = new TFile(filename.c_str(), "RECREATE",
			"CRPropa output data file");
	Ntuple =
			new TNtuple("traj", "CRPropa 3D trajectories",
					"Particle_Type:Initial_Type:Time_Mpc:Position_X_Mpc:Position_Y_Mpc:Position_Z_Mpc:Direction_X:Direction_Y:Direction_Z:Energy_EeV");
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
		Ntuple->Fill(c->current.getId(), c->source.getId(),
				c->getTrajectoryLength() / Mpc,
				c->current.getPosition().getX() / Mpc,
				c->current.getPosition().getY() / Mpc,
				c->current.getPosition().getZ() / Mpc,
				c->current.getDirection().getX(),
				c->current.getDirection().getY(),
				c->current.getDirection().getZ(), c->current.getEnergy() / EeV);
	}
	TThread::UnLock();
}
////////////////////////////////////////////////////////////////////////////////

}// namespace crpropa

#endif // CRPROPA_HAVE_ROOT

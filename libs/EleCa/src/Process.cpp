#include "EleCa/Process.h"

#include <string>
#include <iostream>
#include <cstdlib>

namespace eleca {

void Process::SetName(std::string nm) {
	fname = nm;
}
const std::string &Process::GetName() const {
	return fname;
}

void Process::SetInteractionAngle(double a) {
	fInteractionAngle = a;
}
double Process::GetInteractionAngle() const {
	return fInteractionAngle;
}

void Process::SetLambda(double le) {
	flambda = le;
}
double Process::GetLambda() const {
	return flambda;
}

void Process::SetLimits(double smin, double smax) {
	fsmin = smin;
	fsmax = smax;
}

void Process::SetLimits() {
	SetLimits(fPi, fname);
}

void Process::SetMax(double smax) {
	fsmax = smax;
}
void Process::SetMin(double smin) {
	fsmin = smin;
}
double Process::GetMin() const {
	return fsmin;
}
double Process::GetMax() const {
	return fsmax;
}

void Process::SetCMEnergy(double s) {
	fCMEnergy = s;
}

void Process::SetCMEnergy(Particle p1, Particle pb) {
	fCMEnergy = 2 * p1.GetEnergy() * pb.GetEnergy()
			* (1 - p1.GetBeta() * cos(fInteractionAngle))
			+ p1.GetMass() * p1.GetMass() + pb.GetMass() * pb.GetMass();
}

void Process::SetCMEnergy() {
	fCMEnergy = 2 * fPi.GetEnergy() * fPt.GetEnergy()
			* (1 - fPi.GetBeta() * cos(fInteractionAngle))
			+ fPi.GetMass() * fPi.GetMass() + fPt.GetMass() * fPt.GetMass();

}

double Process::GetCMEnergy() const {
	return fCMEnergy;
}

void Process::SetIncidentParticle(const Particle& p1) {
	fPi = p1;
	SetLimits();
}
void Process::SetTargetParticle(Particle& p1) {
	fPt = p1;
	SetLimits();
}

const Particle &Process::GetIncidentParticle() const {
	return fPi;
}
const Particle &Process::GetTargetParticle() const {
	return fPt;
}

const std::string &Process::GetBackground() const {
	return fback;
}

Process::Process() {
	fname = "";
	SetLimits(0.0, 1.0e23);
	flambda = 0;
	fCMEnergy = 0;
	fInteractionAngle = cPI;
	fback = "ALL";
	fbackdensity = 0;
	feps_inf = eps_ph_inf_global;
	feps_sup = eps_ph_sup_global;
}

Process::Process(Particle& p1, Particle& p2) {
	fPi = p1;
	fPt = p2;
	if (p1.GetType() == 22)
		fname = "PP";
	else if (abs(p1.GetType()) == 11) {
		std::cerr << "NB: by default process set to ICS" << std::endl;
		fname = "ICS";
	} else
		fname = "NONE";
	SetCMEnergy(p1, p2);
	flambda = 0;
	fInteractionAngle = cPI;
	fback = "ALL";
	SetLimits(p1, fname);
	fbackdensity = 0;
	feps_inf = eps_ph_inf_global;
	feps_sup = eps_ph_sup_global;
}

Process::Process(Particle& p1, Particle& p2, std::string name) {
	fname = name;
	SetCMEnergy(p1, p2);
	flambda = 0;
	fInteractionAngle = cPI;
	fPi = p1;
	fPt = p2;
	fback = "ALL";
	SetLimits(p1, fname);
	fbackdensity = 0;
	feps_inf = eps_ph_inf_global;
	feps_sup = eps_ph_sup_global;
}

Process::Process(const Process& proc2) {
	fname = proc2.GetName();
	SetLimits(proc2.GetMin(), proc2.GetMax());
	fCMEnergy = proc2.GetCMEnergy();
	fInteractionAngle = proc2.GetInteractionAngle();
	fPi = proc2.GetIncidentParticle();
	fPt = proc2.GetTargetParticle();
	fback = proc2.GetBackground();
	fbackdensity = 0;
	feps_inf = eps_ph_inf_global;
	feps_sup = eps_ph_sup_global;
}

Process::~Process() {
}

//-----------

void Process::SetBackground(std::string BackRad) {

	fback = BackRad;

	double eps_min = eps_ph_inf_global;
	double eps_max = eps_ph_sup_global;

	if (BackRad == "CMB") {
		eps_min = eps_ph_inf_cmb;
		eps_max = eps_ph_sup_cmb;
	} else if (BackRad == "COB") {
		eps_min = eps_ph_inf_cob;
		eps_max = eps_ph_sup_cob;
	} else if (BackRad == "CIB") {
		eps_min = eps_ph_inf_cib;
		eps_max = eps_ph_sup_cib;
	} else if (BackRad == "CIOB") {
		eps_min = eps_ph_inf_ciob;
		eps_max = eps_ph_sup_ciob;
	} else if (BackRad == "URB") {
		eps_min = eps_ph_inf_urb;
		eps_max = eps_ph_sup_urb;
	}

	feps_inf = eps_min;
	feps_sup = eps_max;

#ifdef DEBUG_ELECA
	std::cout << "eps range set to " << eps_min << " , " << eps_max
	<< std::endl;
#endif

}

void Process::SetLimits(Particle& p1, std::string nameproc) {
	if (p1.GetType() != 22 && p1.GetType() != 11 && p1.GetType() != -11)
		std::cout << "error in type " << p1.GetType() << " != 11 and !=22 "
				<< std::endl;

	if (nameproc == "PP") {
		if (abs(p1.GetType()) == 11)
			std::cout << "\nERROR!! wrong particle or process!! " << nameproc
					<< p1.GetType() << "\n" << std::endl;
		fsmin = 4 * ElectronMass * ElectronMass; //min per theta = 0;
		fsmax = 4 * p1.GetEnergy() * feps_sup * 2.0; //(1.0-cos(fInteractionAngle)); //max per theta=3.14
	}
	if (nameproc == "DPP") {
		if (abs(p1.GetType()) == 11)
			std::cout << "\nERROR!! wrong particle or process!! " << nameproc
					<< p1.GetType() << "\n" << std::endl;
		fsmin = 16 * ElectronMass * ElectronMass;
		fsmax = 4 * p1.GetEnergy() * feps_sup * 2.0;
	}
	if (nameproc == "ICS") {
		if (p1.GetType() == 22 || p1.GetType() == 9)
			std::cout << "\nERROR!! wrong particle or process!! " << nameproc
					<< p1.GetType() << "\n" << std::endl;
		fsmin = p1.GetMass() * p1.GetMass()
				+ 2 * p1.GetEnergy() * feps_inf * (1 - p1.GetBeta());
		fsmax = p1.GetMass() * p1.GetMass()
				+ 2 * p1.GetEnergy() * feps_sup * (1 + p1.GetBeta());
	}
	if (nameproc == "TPP") {
		if (p1.GetType() == 22 || p1.GetType() == 9)
			std::cout << "\nERROR!! wrong particle or process!! " << nameproc
					<< p1.GetType() << "\n" << std::endl;
		fsmin = 10 * ElectronMass * ElectronMass;
		fsmax = 2 * p1.GetEnergy() * feps_sup * (1 + p1.GetBeta())
				+ p1.GetMass() * p1.GetMass();
	}
}

} // namespace

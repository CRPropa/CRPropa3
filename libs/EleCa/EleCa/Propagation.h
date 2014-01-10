#ifndef ELECA_PROPAGATION_H
#define ELECA_PROPAGATION_H

#include "Particle.h"
#include "Process.h"
#include "Common.h"
#include "XLoss_CBR.h"
#include "EnergyLoss.h"
#include "Constants.h"

#include <string>
#include <vector>
#include <fstream>

namespace eleca {

class Propagation {

private:

	double vPPle[1101];
	double vDPPle[1101];
	double vTPPle[1101];
	double vICSle[1101];
	double vEtab[1101];

	std::vector<double> BkgE, BkgA;
	std::string Bkg;
	double fEthr;

public:

	Propagation();

	void SetEthr(double eth) {
		fEthr = eth;
	}

	double GetEthr() {
		return fEthr;
	}

	~Propagation() {
	}

	void WriteOutput(std::ostream &out, Particle &p1,
			std::vector<Particle> &part) const;

	void ReadTables(const std::string &file);
	void InitBkgArray(const std::string &BackRad);

	double GetMeanThetaBFDeflection(double Bin, double Ein, int ptype,
			double Lin) const;
	double GetLambdaTab(const Process &proc, const std::string &procName) const;
	double ExtractMinDist(Process &proc, int type, double R, double R2,
			std::vector<double> &Etarget) const;
	std::vector<double> GetEtarget(Process &proc,
			const Particle &particle) const;
	void Propagate(Particle &curr_particle,
			std::vector<Particle> &ParticleAtMatrix,
			std::vector<Particle> &ParticleAtGround) const;
	double ExtractPhotonEnergyMC(double z, Process &proc) const;
	double ShootPhotonEnergyMC(double z) const;
	double ShootPhotonEnergyMC(double Emin, double z) const;
	void SetInitVar(std::vector<std::vector<double> > bk,
			std::vector<std::vector<double> > *le) const;

};

Propagation::Propagation() {
	fEthr = 1e16;
}

void Propagation::ReadTables(const std::string &filename) {

#ifdef DEBUG_ELECA
	std::cout << filename << std::endl;
#endif

	std::ifstream fin(filename.c_str());

	if (!fin.is_open()) {
		std::cerr << "Unable to open lambda_table file: " << filename
				<< " ! exiting... ";
		return;
	}

	int k = 0;
	double Etab, PPle, ICSle, DPPle, TPPle;
	while (fin.good()) {
		fin >> Etab >> PPle >> ICSle >> DPPle >> TPPle;
		vEtab[k] = Etab;
		vPPle[k] = PPle;
		vICSle[k] = ICSle;
		vDPPle[k] = DPPle;
		vTPPle[k] = TPPle;
		k++;
	}

	if (k != 1101)
		std::cerr << "Failed to read lambda_table file: " << filename
				<< "! only " << k << " entries, expected 1101!";
}

void Propagation::InitBkgArray(const std::string &BackRad) {
	// Routine to build the array of cumulative distribution of
	// background photons

	Bkg = BackRad;
	BkgE.resize(POINTS_VERY_FEW);
	BkgA.resize(POINTS_VERY_FEW);

	if (BackRad == "CMB") {
		double de = pow((double) eps_ph_sup_cmb / eps_ph_inf_cmb,
				1. / POINTS_VERY_FEW);
		double e = eps_ph_inf_cmb;
		for (size_t i = 0; i < POINTS_VERY_FEW; i++) {
			BkgE[i] = e;
			BkgA[i] = CMBR(e);
			e *= de;
		}
	}

	if (BackRad == "CIOB") {
		double de = pow((double) eps_ph_sup_ciob / eps_ph_inf_ciob,
				1. / POINTS_VERY_FEW);
		double e = eps_ph_inf_ciob;
		for (size_t i = 0; i < POINTS_VERY_FEW; i++) {
			BkgE[i] = e;
			BkgA[i] = CIOBR(e);
			e *= de;
		}
	}

	if (BackRad == "URB") {
		double de = pow((double) eps_ph_sup_urb / eps_ph_inf_urb,
				1. / POINTS_VERY_FEW);
		double e = eps_ph_inf_urb;
		for (size_t i = 0; i < POINTS_VERY_FEW; i++) {
			BkgE[i] = e;
			BkgA[i] = URB(e);
			e *= de;
		}
	}

	if (BackRad == "URB") {
		double de = pow((double) eps_ph_sup_global / eps_ph_inf_global,
				(double) 1. / POINTS_VERY_FEW);
		double e = eps_ph_inf_global;
		for (size_t i = 0; i < POINTS_VERY_FEW; i++) {
			BkgE[i] = e;
			BkgA[i] = CBR(e);
			e *= de;
		}
	}

	// cumulate
	for (size_t i = 1; i < POINTS_VERY_FEW; i++) {
		BkgA[i] += BkgA[i - 1];
	}

	// normalize
	double a = 1.0 / BkgA[POINTS_VERY_FEW - 1];
	for (size_t i = 0; i < POINTS_VERY_FEW; i++) {
		BkgA[i] *= a;
#ifdef DEBUG_ELECA
		std::cout << BkgE[i] << " - " << BkgA[i] << std::endl;
#endif
	}
}

double Propagation::GetMeanThetaBFDeflection(double Bin, double Ein, int ptype,
		double Lin) const {
	//from D. Hooper, S. Sarkar, M. Taylor, arXiv: 0608085, 2006

	if (Bin == 0 || Ein == 0)
		return 0;
	if (ptype == 22)
		return 0;

	double lcoher = 1;

	return 0.8 * (1.0e20 / Ein) * sqrt(Lin / 10 * lcoher) * (Bin / 1.0e-9)
			/ 180.0 * 3.1415927;
}

double Propagation::ExtractMinDist(Process &proc, int type, double R, double R2,
		std::vector<double> &Etarget) const {

	double min_dist1 = 0;
	double min_dist2 = 0;
	Process proc1(proc);
	Process proc2(proc);
	double tmp_lambda1 = 0;
	double tmp_lambda2 = 0;
	Particle pt;
	pt.SetType(0);
	pt.Setz(proc.GetIncidentParticle().Getz());

	if (type == 22) {

		proc1.SetName("PP");
		pt.SetEnergy(Etarget[0]);
		proc1.SetTargetParticle(pt);
		proc1.SetCMEnergy();

		tmp_lambda1 = GetLambdaTab(proc1, "PP");

		min_dist1 = -tmp_lambda1 * log(R);

		proc2.SetName("DPP");
		pt.SetEnergy(Etarget[1]);
		proc2.SetTargetParticle(pt);
		proc2.SetCMEnergy();
		tmp_lambda2 = GetLambdaTab(proc2, "DPP");
		min_dist2 = -tmp_lambda2 * log(R2);

#ifdef DEBUG_ELECA
		std::cerr << "comparing 2 mindists: " << min_dist1 << "("
		<< tmp_lambda1 << ") vs " << min_dist2 << " ( "
		<< tmp_lambda2 << ") " << std::endl;
#endif

		if (min_dist2 < min_dist1) {
			min_dist1 = min_dist2;
			proc.SetName("DPP");
			pt.SetEnergy(Etarget[1]);
			proc.SetTargetParticle(pt);
			proc.SetCMEnergy();
		} else {
			proc.SetName("PP");
			pt.SetEnergy(Etarget[0]);
			proc.SetTargetParticle(pt);
			proc.SetCMEnergy();
		}
	}    //end if type 0
	else if (abs(type) == 11) {

		proc1.SetName("ICS");
		pt.SetEnergy(Etarget[0]);
		proc1.SetTargetParticle(pt);
		tmp_lambda1 = GetLambdaTab(proc1, "ICS");
		min_dist1 = -tmp_lambda1 * log(R);

		proc2.SetName("TPP");
		pt.SetEnergy(Etarget[1]);
		proc2.SetTargetParticle(pt);
		tmp_lambda2 = GetLambdaTab(proc2, "TPP");
		min_dist2 = -tmp_lambda2 * log(R2);

#ifdef DEBUG_ELECA
		std::cerr << "comparing 2 mindists: " << min_dist1 << "("
		<< tmp_lambda1 << ") vs " << min_dist2 << " ( "
		<< tmp_lambda2 << ") " << std::endl;
#endif

		if (min_dist2 < min_dist1) {
			min_dist1 = min_dist2;
			proc.SetName("TPP");
			pt.SetEnergy(Etarget[1]);
			proc.SetTargetParticle(pt);
			proc.SetCMEnergy();
		} else {
			proc.SetName("ICS");
			pt.SetEnergy(Etarget[0]);
			proc.SetTargetParticle(pt);
			proc.SetCMEnergy();
		}
	}    //else e+/e-
	else
		std::cerr << "something wrong in particle type ( " << type
				<< ". Propagation of photons and e+/e- is the only allowed.)"
				<< std::endl;

	return min_dist1;
}

double Propagation::GetLambdaTab(const Process &proc,
		const std::string &procName) const {

	double E1 = proc.GetIncidentParticle().GetEnergy();
	double z = proc.GetIncidentParticle().Getz();
	double res = 0;

	double E0taborg = vEtab[0];

	double dEtab = log10(vEtab[0]) - log10(vEtab[1]);
	double evolution = GetEvolution(proc.GetTargetParticle().GetEnergy(), z);
	int i = (int) ((log10(E0taborg) - log10(E1 * (1 + z))) / dEtab);

	if (i < 0) {
		std::cout << "WARNING!! GetLambdaTab in " << procName << " : i= " << i
				<< " <0! E1*(1+z) =   " << E1 << "* (1 + " << z << ") < "
				<< E0taborg << ".. returning lambda[0];" << std::endl;
	}

	else if (i >= 1001) {
		std::cout << "WARNING!! GetLambdaTab in " << procName << " : i>= "
				<< 1001 << " ! E1*(1+z) =   " << E1 << "* (1 + " << z
				<< ") .. returning lambda[nentries];" << std::endl;

	} else {
		if (procName == "PP")
			res = vPPle[i];
		else if (procName == "DPP")
			res = vDPPle[i];
		else if (procName == "ICS")
			res = vICSle[i];
		else if (procName == "TPP")
			res = vTPPle[i];
	}

	if (evolution != 0) {
		if (res / evolution < 0)
			std::cerr
					<< "ERROR UNPHYSICAL SOLUTION!! CHECK HERE LAMBDA OR EVOLUTION!!"
					<< std::endl;
		return res / evolution;
	}
	std::cerr << "warning!! evolution ==0 " << std::endl;
	return 0;
}

double Propagation::ShootPhotonEnergyMC(double z) const {
	// Routine for the MC sampling of background photon energy

//	double interpolate(double x, const std::vector<double> &X,
//			const std::vector<double> &Y) {
//		std::vector<double>::const_iterator it = std::upper_bound(X.begin(),
//				X.end(), x);
//		if (it == X.begin())
//			return Y.front();
//		if (it == X.end())
//			return Y.back();
//
//		size_t i = it - X.begin() - 1;
//		return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
//	}

	double h = Uniform(0, 1);
	for (int i = 0; i < POINTS_VERY_FEW; i++) {
		if (h < BkgA[i]) {
			return BkgE[i] * (1. + z);
			break;
		}
	}
#ifdef DEBUG_ELECA
	std::cout << "ShootPhotonEnergyMC. z = " << z << " h: " << h << " => 0"
	<< std::endl;
#endif

	return 0.;
}

double Propagation::ShootPhotonEnergyMC(double Emin, double z) const {
	// Routine for the MC sampling of background photon energy
	std::vector<double>::const_iterator it;

	// find lowest energy bin

	it = std::lower_bound(BkgE.begin(), BkgE.end(), Emin);
	size_t iE;
	if (it == BkgE.begin())
		iE = 0;
	else if (it == BkgE.end())
		iE = BkgE.size() - 1;
	else
		iE = it - BkgE.begin();

	// random number in selected range
	double h = Uniform(BkgA[iE], 1);

	// find energy for random number
	it = std::upper_bound(BkgA.begin(), BkgA.end(), h);
	if (it == BkgA.begin())
		return BkgE.front();
	else if (it == BkgA.end())
		return BkgE.back();
	else
		return BkgE[it - BkgA.begin()];
}

std::vector<double> Propagation::GetEtarget(Process &proc,
		const Particle &particle) const {

	std::vector<double> Etarget;
	double Etarget_tmp = 0;
	double smintmp = 0;
	double z_curr = particle.Getz();
	double Energy = particle.GetEnergy();
	int pType = particle.GetType();
	if (pType == 22) {
		proc.SetName("PP");
		proc.SetLimits();
		smintmp = proc.GetMin();
		Etarget_tmp = ShootPhotonEnergyMC(ElectronMass * ElectronMass / Energy,
				z_curr);
		Etarget.push_back(Etarget_tmp);

		proc.SetName("DPP");
		proc.SetLimits();
		smintmp = proc.GetMin();
		Etarget_tmp = ShootPhotonEnergyMC(smintmp / (4.0 * Energy), z_curr);
		Etarget.push_back(Etarget_tmp);
	}

	else if (abs(pType) == 11) {
		proc.SetName("ICS");
		proc.SetLimits();
		smintmp = proc.GetMin();
		Etarget_tmp = ShootPhotonEnergyMC(smintmp / (4.0 * Energy), z_curr);
		Etarget.push_back(Etarget_tmp);

		proc.SetName("TPP");
		proc.SetLimits();
		smintmp = proc.GetMin();
		Etarget_tmp = ShootPhotonEnergyMC(smintmp / (4.0 * Energy), z_curr);
		Etarget.push_back(Etarget_tmp);
	}    //end e/e
	else
		std::cerr << "something wrong in particle type ( " << pType
				<< ". Propagation of photons and e+/e- is the only allowed.)"
				<< std::endl;

	if (Etarget.size() != 2) {
		std::cout << "something wrong with the Etarget!! " << std::endl;
		exit(0);
	}

	return Etarget;
}

double Propagation::ExtractPhotonEnergyMC(double z, Process &proc) const {
	double esoft = 0;
//double snew = 0;
	double emin = proc.GetMin();
	Particle pi = proc.GetIncidentParticle();
	Particle pb = proc.GetTargetParticle();

	double Epi = pi.GetEnergy();
	double m = pi.GetMass();
	esoft = ShootPhotonEnergyMC(emin / (4.0 * Epi), z);
	//snew = 4 * Epi * esoft + m * m;
	pb.SetEnergy(esoft);
	proc.SetTargetParticle(pb);
	proc.SetCMEnergy();
	return esoft;
}

void Propagation::WriteOutput(std::ostream &out, Particle &p1,
		std::vector<Particle> &ParticleAtGround) const {
	double Bfield = 0;
	size_t NsecG = ParticleAtGround.size();

	out << fEthr << " " << Bfield / 1e-9 << " " << p1.GetEnergy() << " "
			<< p1.Getz() << " " << NsecG;
	for (int i = 0; i < NsecG; ++i) {
		Particle &p = ParticleAtGround[i];
		out << "  " << p.GetWeigth() << " " << p.GetEnergy() << " "
				<< p.GetType();
	}
	out << std::endl;
}
//
//void Propagation::Spectrum(std::vector<double> &spectrum) const {
//	double emin = 7.0;
//	double dE = (24.0 - 7.0) / 170.0;
//	size_t ipos = 0;
//	size_t NsecG = ParticleAtGround.size();
//
//	for (int h = 0; h < NsecG; ++h) {
//		ipos = (int) ((log10(EGround.at(h)) - emin) / dE);
//		if (typeGround.at(h) == 22)
//			fdN[ipos] += wGround.at(h);
//	}
//}
//
//void Propagation::AddSpectrum(std::vector<double> &spectrum) const {
//	double emin = 7.0;
//	double dE = (24.0 - 7.0) / 170.0;
//	size_t ipos = 0;
//	size_t NsecG = ParticleAtGround.size();
//
//	for (int h = 0; h < NsecG; ++h) {
//		ipos = (int) ((log10(EGround.at(h)) - emin) / dE);
//		if (typeGround.at(h) == 22)
//			fdN[ipos] += wGround.at(h);
//	}
//}

void Propagation::Propagate(Particle &curr_particle,
		std::vector<Particle> &ParticleAtMatrix,
		std::vector<Particle> &ParticleAtGround) const {

	double theta_deflBF = 0.0;
	double BNorm = curr_particle.GetB();

	double zin = curr_particle.Getz();
	double Ein = curr_particle.GetEnergy();
	int type = curr_particle.Getz();

	int wi_last = curr_particle.GetWeigth();

	double z_curr = zin;
	double Ecurr = Ein;

	bool interacted = 0;
	double min_dist = 1e12;
	double walkdone = 0;

	double E1 = 0;
	double E2 = 0;
	double E3 = 0;

	double stepsize = 0;
	double Elast = 0;

	double R = Uniform(0.0, 1.0);
	double R2 = Uniform(0.0, 1.0);
	bool fast = 0;

	Process proc;
	proc.SetIncidentParticle(curr_particle);
	proc.SetBackground(Bkg);

	std::vector<double> EtargetAll = GetEtarget(proc, curr_particle);
#ifdef DEBUG_ELECA
	std::cout << "GetEtarget: " << EtargetAll[0] << " " << EtargetAll[1] << std::endl;
#endif

	min_dist = ExtractMinDist(proc, curr_particle.GetType(), R, R2, EtargetAll);
#ifdef DEBUG_ELECA
	std::cout << "ExtractMinDist " << min_dist << std::endl;
#endif

	interacted = 0;
	double dz = 0;
	double zpos = zin;

	double corrB_factor = 0;
	double realpath = 0;

	double min_dist_last = min_dist;
#ifdef DEBUG_ELECA
	std::cout << "starting propagation... min_dist_last: " << min_dist_last << std::endl;
#endif

	while (!interacted) {

		proc.SetInteractionAngle(cPI);
		theta_deflBF = 0;
		realpath = 0.1 * min_dist;

		theta_deflBF = GetMeanThetaBFDeflection(BNorm,
				curr_particle.GetEnergy(), curr_particle.GetType(), min_dist);
		corrB_factor = cos(theta_deflBF);

		stepsize = realpath * corrB_factor;
		dz = Mpc2z(stepsize);

#ifdef DEBUG_ELECA
		if (BNorm > 0)
		std::cout << "z_curr " << z_curr << ", type: "
		<< curr_particle.GetType() << ", Ecurr "
		<< curr_particle.GetEnergy() << "  = " << min_dist
		<< " Mpc -->  theta: " << theta_deflBF * 180.0 / 3.1415927
		<< "  " << corrB_factor << std::endl;

		std::cout << "done : " << walkdone << " + " << realpath << " = "
		<< walkdone + realpath << " has to be:  " << min_dist
		<< " Mpc. Now @ z= " << zpos << " vs nexz = "
		<< zpos - Mpc2z(stepsize) << std::endl;
#endif

		if ((walkdone + realpath) > min_dist) {
#ifdef DEBUG_ELECA
			std::cout
			<< " walkdone + realpath > min_dist: correcting realpath from"
			<< realpath << " to " << min_dist - walkdone
			<< " and the stepsize is changed from " << dz << " to ";
			realpath = min_dist - walkdone;
			stepsize = realpath * corrB_factor;
			dz = Mpc2z(stepsize);
			std::cout << dz << std::endl;
#endif
			interacted = 1;
		}

		if (zpos - dz <= 0) {
			dz = zpos;
			stepsize = z2Mpc(dz);
			realpath = stepsize / corrB_factor;
		}

		zpos -= dz;
		walkdone += realpath;
		Elast = Ecurr;

		double adiab = Ecurr
				- AdiabaticELoss(zpos + Mpc2z(realpath), zpos + dz, Ecurr);

		if (type == 0 || type == 22)
			Ecurr = EnergyLoss1D(Ecurr, zpos + Mpc2z(realpath), zpos, 0);
		else
			Ecurr = EnergyLoss1D(Ecurr, zpos + Mpc2z(realpath), zpos, BNorm);

		z_curr = zpos;

		curr_particle.Setz(z_curr);
		curr_particle.SetEnergy(Ecurr);

		if (z_curr > 0 && Ecurr < fEthr) {
			return;
		}
		if (z_curr <= 0) {
			ParticleAtGround.push_back(curr_particle);
			return;
		}

		proc.SetIncidentParticle(curr_particle);
		proc.SetCMEnergy();
		proc.SetLimits();
		//      std::vector<double> EtargetAll=GetEtarget(proc,curr_particle);
		min_dist = ExtractMinDist(proc, curr_particle.GetType(), R, R2,
				EtargetAll);
	} //end while

	if (interacted == 1) {
#ifdef DEBUG_ELECA
		std::cerr << "******** producing secondary particles according to "
		<< proc.GetName() << " process *********** " << std::endl;
#endif
		if (proc.GetName() == "PP") {

			E1 = ExtractPPSecondariesEnergy(proc);

			if (E1 == 0 || E1 == Ecurr)
				std::cerr << "ERROR in PP process:  E : " << Ecurr << "  " << E1
						<< " " << std::endl;

			if (E1 > fEthr) {
				Particle pp(11, E1, z_curr);
				pp.SetWeigth(wi_last);
				ParticleAtMatrix.push_back(pp);
			}
			if (Ecurr - E1 > fEthr) {
				Particle pe(-11, Ecurr - E1, z_curr);
				pe.SetWeigth(wi_last);
				ParticleAtMatrix.push_back(pe);
			}

			return;
		} //if PP
		else if (proc.GetName() == "DPP") {
			E1 = (Ecurr - 2 * ElectronMass) / 2.0;
			if (E1 == 0)
				std::cerr << "ERROR in DPP process E : " << E1 << std::endl;

			if (E1 > fEthr) {
				Particle pp(11, E1, z_curr);
				if (fast == 1)
					pp.SetWeigth(wi_last * 2);
				else {
					pp.SetWeigth(wi_last);
					Particle pe(-11, E1, z_curr);
					pe.SetWeigth(wi_last);
					ParticleAtMatrix.push_back(pe);
				}
				ParticleAtMatrix.push_back(pp);
			}

			return;
		} //end if DPP
		else if (proc.GetName() == "ICS") {

			E1 = ExtractICSSecondariesEnergy(proc);
			E2 = Ecurr - E1;
			if (E1 == 0 || E2 == 0)
				std::cerr << "ERROR in ICS process E : " << E1 << " " << E2
						<< std::endl;

			if (E1 > fEthr) {
				Particle pp(curr_particle.GetType(), E1, z_curr);
				pp.SetWeigth(wi_last);
				ParticleAtMatrix.push_back(pp);
			}
			if (E2 > fEthr) {
				Particle pg(22, E2, z_curr);
				pg.SetWeigth(wi_last);
				ParticleAtMatrix.push_back(pg);
			}

			return;
		} //if ics
		else if (proc.GetName() == "TPP") {
			E1 = E2 = ExtractTPPSecondariesEnergy(proc);
			E3 = Ecurr - E1 - E2;
			if (E1 == 0 || E2 == 0 || E3 == 0)
				std::cerr << "ERROR in TPP process E : " << E1 << " " << E2
						<< std::endl;

			if (E1 > fEthr) {
				Particle pp(11, E1, z_curr);
				if (fast == 1)
					pp.SetWeigth(wi_last * 2);
				else {
					pp.SetWeigth(wi_last);
					Particle pe(-11, E1, z_curr);
					pe.SetWeigth(wi_last);
					ParticleAtMatrix.push_back(pe);
				}
				ParticleAtMatrix.push_back(pp);
			}
			if (E3 > fEthr) {
				Particle psc(curr_particle.GetType(), E3, z_curr);
				psc.SetWeigth(wi_last);
				ParticleAtMatrix.push_back(psc);
			}
			return;
		}
	}

	return;

}

} // namespace eleca

#endif // ELECA_PROPAGATION_H

#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"

#include "kiss/convert.h"
#include "kiss/logger.h"
#include "sophia.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {

PhotoPionProduction::PhotoPionProduction(ref_ptr<PhotonField> field, bool photons, bool neutrinos, bool electrons, bool antiNucleons, double l, bool redshift) {
	havePhotons = photons;
	haveNeutrinos = neutrinos;
	haveElectrons = electrons;
	haveAntiNucleons = antiNucleons;
	haveRedshiftDependence = redshift;
	limit = l;
	setPhotonField(field);
}

void PhotoPionProduction::setPhotonField(ref_ptr<PhotonField> field) {
	photonField = field;
	std::string fname = photonField->getFieldName();
	if (haveRedshiftDependence) {
		if (photonField->hasRedshiftDependence() == false){
			std::cout << "PhotoPionProduction: tabulated redshift dependence not needed for " + fname + ", switching off" << std::endl;
			haveRedshiftDependence = false;
		}
		else {
			KISS_LOG_WARNING << "PhotoPionProduction: You are using the 2-dimensional tabulated redshift evolution, which is not available for other interactions. To be consistent across all interactions you may deactivate this <setHaveRedshiftDependence(False)>.";
		}
	}
	
	setDescription("PhotoPionProduction: " + fname);
	if (haveRedshiftDependence){
		initRate(getDataPath("PhotoPionProduction/rate_" + fname.replace(0, 3, "IRBz") + ".txt"));
	}
	else
		initRate(getDataPath("PhotoPionProduction/rate_" + fname + ".txt"));
}

void PhotoPionProduction::setHavePhotons(bool b) {
	havePhotons = b;
}
	
void PhotoPionProduction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void PhotoPionProduction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void PhotoPionProduction::setHaveAntiNucleons(bool b) {
	haveAntiNucleons = b;
}

void PhotoPionProduction::setHaveRedshiftDependence(bool b) {
	haveRedshiftDependence = b;
	setPhotonField(photonField);
}

void PhotoPionProduction::setLimit(double l) {
	limit = l;
}

void PhotoPionProduction::initRate(std::string filename) {
	// clear previously loaded tables
	tabLorentz.clear();
	tabRedshifts.clear();
	tabProtonRate.clear();
	tabNeutronRate.clear();

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("PhotoPionProduction: could not open file " + filename);

	if (haveRedshiftDependence) {
		double zOld = -1, aOld = -1;
		while (infile.good()) {
			if (infile.peek() == '#') {
				infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				continue;
			}
			double z, a, b, c;
			infile >> z >> a >> b >> c;
			if (!infile)
				break;
			if (z > zOld) {
				tabRedshifts.push_back(z);
				zOld = z;
			}
			if (a > aOld) {
				tabLorentz.push_back(pow(10, a));
				aOld = a;
			}
			tabProtonRate.push_back(b / Mpc);
			tabNeutronRate.push_back(c / Mpc);
		}
	} else {
		while (infile.good()) {
			if (infile.peek() == '#') {
				infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				continue;
			}
			double a, b, c;
			infile >> a >> b >> c;
			if (!infile)
				break;
			tabLorentz.push_back(pow(10, a));
			tabProtonRate.push_back(b / Mpc);
			tabNeutronRate.push_back(c / Mpc);
		}
	}

	infile.close();
}

double PhotoPionProduction::nucleonMFP(double gamma, double z, bool onProton) const {
	const std::vector<double> &tabRate = (onProton)? tabProtonRate : tabNeutronRate;

	// scale nucleus energy instead of background photon energy
	gamma *= (1 + z);
	if (gamma < tabLorentz.front() or (gamma > tabLorentz.back()))
		return std::numeric_limits<double>::max();

	double rate;
	if (haveRedshiftDependence)
		rate = interpolate2d(z, gamma, tabRedshifts, tabLorentz, tabRate);
	else
		rate = interpolate(gamma, tabLorentz, tabRate) * photonField->getRedshiftScaling(z);

	// cosmological scaling
	rate *= pow_integer<2>(1 + z);

	return 1. / rate;
}

double PhotoPionProduction::nucleiModification(int A, int X) const {
	if (A == 1)
		return 1.;
	if (A <= 8)
		return 0.85 * pow(X, 2. / 3.);
	return 0.85 * X;
}

void PhotoPionProduction::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	double z = candidate->getRedshift();
	// the loop is processed at least once for limiting the next step
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (!isNucleus(id))
			return;

		// find interaction with minimum random distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();
		double meanFreePath;
		double totalRate = 0;
		bool onProton = true; // interacting particle: proton or neutron

		int A = massNumber(id);
		int Z = chargeNumber(id);
		int N = A - Z;
		double gamma = candidate->current.getLorentzFactor();

		// check for interaction on protons
		if (Z > 0) {
			meanFreePath = nucleonMFP(gamma, z, true) / nucleiModification(A, Z);
			randDistance = -log(random.rand()) * meanFreePath;
			totalRate += 1. / meanFreePath;
		}
		// check for interaction on neutrons
		if (N > 0) {
			meanFreePath = nucleonMFP(gamma, z, false) / nucleiModification(A, N);
			totalRate += 1. / meanFreePath;
			double d = -log(random.rand()) * meanFreePath;
			if (d < randDistance) {
				randDistance = d;
				onProton = false;
			}
		}

		// check if interaction does not happen
		if (step < randDistance) {
			if (totalRate > 0.)
				candidate->limitNextStep(limit / totalRate);
			return;
		}

		// interact and repeat with remaining step
		performInteraction(candidate, onProton);
		step -= randDistance;
	} while (step > 0);
}

void PhotoPionProduction::performInteraction(Candidate *candidate, bool onProton) const {
	int id = candidate->current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);
	double E = candidate->current.getEnergy();
	double EpA = E / A;
	double z = candidate->getRedshift();

	// SOPHIA simulates interactions only for protons / neutrons.
	// For anti-protons / neutrons assume charge symmetry and change all
	// interaction products from particle <--> anti-particle (sign)
	int sign = (id > 0) ? 1 : -1;

	// check if below SOPHIA's energy threshold
	double E_threshold = (photonField->getFieldName() == "CMB") ? 3.72e18 * eV : 5.83e15 * eV;
	if (EpA * (1 + z) < E_threshold)
		return;

	// SOPHIA - input:
	int nature = 1 - static_cast<int>(onProton);  // 0=proton, 1=neutron
	double Ein = EpA / GeV;  // GeV is the SOPHIA standard unit
	double eps = sampleEps(onProton, EpA, z) / GeV;  // GeV for SOPHIA

	// SOPHIA - output:
	double outputEnergy[5][2000];  // [GeV/c, GeV/c, GeV/c, GeV, GeV/c^2]
	int outPartID[2000];
	int nParticles;

#pragma omp critical
	{
		sophiaevent_(nature, Ein, eps, outputEnergy, outPartID, nParticles);
	}

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	std::vector<int> pnType;  // filled with either 13 (proton) or 14 (neutron)
	std::vector<double> pnEnergy;  // corresponding energies of proton or neutron
	if (nParticles == 0)
		return;
	for (int i = 0; i < nParticles; i++) { // loop over out-going particles
		double Eout = outputEnergy[3][i] * GeV; // only the energy is used; could be changed for more detail
		int pType = outPartID[i];
		switch (pType) {
		case 13: // proton
		case 14: // neutron
			// proton and neutron data is taken to determine primary particle in a later step
			pnType.push_back(pType);
			pnEnergy.push_back(Eout);
			break;
		case -13: // anti-proton
		case -14: // anti-neutron
			if (haveAntiNucleons)
				try
				{
					candidate->addSecondary(-sign * nucleusId(1, 14 + pType), Eout, pos, 1., interactionTag);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (anti-nucleon production)\n" << "Something went wrong in the PhotoPionProduction\n"<< "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			break;
		case 1: // photon
			if (havePhotons)
				candidate->addSecondary(22, Eout, pos, 1., interactionTag);
			break;
		case 2: // positron
			if (haveElectrons)
				candidate->addSecondary(sign * -11, Eout, pos, 1., interactionTag);
			break;
		case 3: // electron
			if (haveElectrons)
				candidate->addSecondary(sign * 11, Eout, pos, 1., interactionTag);
			break;
		case 15: // nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * 12, Eout, pos, 1., interactionTag);
			break;
		case 16: // anti-nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * -12, Eout, pos, 1., interactionTag);
			break;
		case 17: // nu_mu
			if (haveNeutrinos)
				candidate->addSecondary(sign * 14, Eout, pos, 1., interactionTag);
			break;
		case 18: // anti-nu_mu
			if (haveNeutrinos)
				candidate->addSecondary(sign * -14, Eout, pos, 1., interactionTag);
			break;
		default:
			throw std::runtime_error("PhotoPionProduction: unexpected particle " + kiss::str(pType));
		}
	}
	double maxEnergy = *std::max_element(pnEnergy.begin(), pnEnergy.end());  // criterion for being declared primary
	for (int i = 0; i < pnEnergy.size(); ++i) {
		if (pnEnergy[i] == maxEnergy) {  // nucleon is primary particle
			if (A == 1) {
				// single interacting nucleon
				candidate->current.setEnergy(pnEnergy[i]);
				try
				{
					candidate->current.setId(sign * nucleusId(1, 14 - pnType[i]));
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (primary particle, A==1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			} else {
				// interacting nucleon is part of nucleus: it is emitted from the nucleus
				candidate->current.setEnergy(E - EpA);
				try
				{
					candidate->current.setId(sign * nucleusId(A - 1, Z - int(onProton)));
					candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos, 1., interactionTag);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (primary particle, A!=1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			}
		} else {  // nucleon is secondary proton or neutron
			candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos, 1., interactionTag);
		}
	}
}

double PhotoPionProduction::lossLength(int id, double gamma, double z) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;

	double lossRate = 0;
	if (Z > 0)
		lossRate += 1 / nucleonMFP(gamma, z, true) * nucleiModification(A, Z);
	if (N > 0)
		lossRate += 1 / nucleonMFP(gamma, z, false) * nucleiModification(A, N);

	// approximate the relative energy loss
	// - nucleons keep the fraction of mass to delta-resonance mass
	// - nuclei lose the energy 1/A the interacting nucleon is carrying
	double relativeEnergyLoss = (A == 1) ? 1 - 938. / 1232. : 1. / A;
	lossRate *= relativeEnergyLoss;

	// scaling factor: interaction rate --> energy loss rate
	lossRate *= (1 + z);

	return 1. / lossRate;
}

SophiaEventOutput PhotoPionProduction::sophiaEvent(bool onProton, double Ein, double eps) const {
	// SOPHIA - input:
	int nature = 1 - static_cast<int>(onProton);  // 0=proton, 1=neutron
	Ein /= GeV;  // GeV is the SOPHIA standard unit
	eps /= GeV;  // GeV for SOPHIA

	// SOPHIA - output:
	double outputEnergy[5][2000];  // [Px GeV/c, Py GeV/c, Pz GeV/c, E GeV, m0 GeV/c^2]
	int outPartID[2000];
	int nParticles;

	sophiaevent_(nature, Ein, eps, outputEnergy, outPartID, nParticles);

	// convert SOPHIA IDs to PDG naming convention & create particles
	SophiaEventOutput output;
	output.nParticles = nParticles;
	for (int i = 0; i < nParticles; ++i) {
		int id = 0;
		int partType = outPartID[i];
		switch (partType) {
			case 13:  // proton
			case 14:  // neutron
				id = nucleusId(1, 14 - partType);
				break;
			case -13:  // anti-proton
			case -14:  // anti-neutron
				id = -nucleusId(1, 14 + partType);
				break;
			case 1:  // photon
				id = 22;
				break;
			case 2:  // positron
				id = -11;
				break;
			case 3:  // electron
				id = 11;
				break;
			case 15:  // nu_e
				id = 12;
				break;
			case 16:  // anti-nu_e
				id = -12;
				break;
			case 17:  // nu_mu
				id = 14;
				break;
			case 18:  // anti-nu_mu
				id = -14;
				break;
			default:
				throw std::runtime_error("PhotoPionProduction: unexpected particle " + kiss::str(partType));
		}
		output.energy.push_back(outputEnergy[3][i] * GeV); // only the energy is used; could be changed for more detail
		output.id.push_back(id);
	}
	return output;
}

double PhotoPionProduction::sampleEps(bool onProton, double E, double z) const {
	// sample eps between epsMin ... epsMax
	double Ein = E / GeV;
	double epsMin = std::max(photonField -> getMinimumPhotonEnergy(z) / eV, epsMinInteraction(onProton, Ein));
	double epsMax = photonField -> getMaximumPhotonEnergy(z) / eV;
	double pEpsMax = probEpsMax(onProton, Ein, z, epsMin, epsMax);

	Random &random = Random::instance();
	for (int i = 0; i < 1000000; i++) {
		double eps = epsMin + random.rand() * (epsMax - epsMin);
		double pEps = probEps(eps, onProton, Ein, z);
		if (random.rand() * pEpsMax < pEps)
			return eps * eV;
	}
	throw std::runtime_error("error: no photon found in sampleEps, please make sure that photon field provides photons for the interaction by adapting the energy range of the tabulated photon field.");
}

double PhotoPionProduction::epsMinInteraction(bool onProton, double Ein) const {
	// labframe energy of least energetic photon where PPP can occur
	// this kind-of ties samplingEps to the PPP and SOPHIA
	const double m = mass(onProton);
	const double p = momentum(onProton, Ein);
	double epsMin = 1.e9 * (1.1646 - m * m) / 2. / (Ein + p); // eV
	return epsMin;
}

double PhotoPionProduction::probEpsMax(bool onProton, double Ein, double z, double epsMin, double epsMax) const {
	// find pEpsMax by testing photon energies (eps) for their interaction
	// probabilities (p) in order to find the maximum (max) probability
	const int nrSteps = 100;
	double pEpsMaxTested = 0.;
	double step = 0.;
	if (sampleLog){
		// sample in logspace with stepsize that is at max Δlog(E/eV) = 0.01 or otherwise dep. on size of energy range with nrSteps+1 steps log. equidis. spaced
		step = std::min(0.01, std::log10(epsMax / epsMin) / nrSteps);
	} else
		step = (epsMax - epsMin) / nrSteps;

	double epsDummy = 0.;
	int i = 0;
	while (epsDummy < epsMax) {
		if (sampleLog)
			epsDummy = epsMin * pow(10, step * i);
		else
			epsDummy = epsMin + step * i;
		double p = probEps(epsDummy, onProton, Ein, z);
		if(p > pEpsMaxTested)
			pEpsMaxTested = p;
		i++;
	}
	// the following factor corrects for only trying to find the maximum on nrIteration photon energies
	// the factor should be determined in convergence tests
	double pEpsMax = pEpsMaxTested * correctionFactor;
	return pEpsMax;
}

double PhotoPionProduction::probEps(double eps, bool onProton, double Ein, double z) const {
	// probEps returns "probability to encounter a photon of energy eps", given a primary nucleon
	// note, probEps does not return a normalized probability [0,...,1]
	double photonDensity = photonField->getPhotonDensity(eps * eV, z) * ccm / eps;
	if (photonDensity != 0.) {
		const double p = momentum(onProton, Ein);
		const double sMax = mass(onProton) * mass(onProton) + 2. * eps * (Ein + p) / 1.e9;
		if (sMax <= sMin())
			return 0;
		double sIntegr = gaussInt([this, onProton](double s) { return this->functs(s, onProton); }, sMin(), sMax);
		return photonDensity * sIntegr / eps / eps / p / 8. * 1.e18 * 1.e6;
	}
	return 0;
}

double PhotoPionProduction::momentum(bool onProton, double Ein) const {
	const double m = mass(onProton);
	const double momentumHadron = sqrt(Ein * Ein - m * m);  // GeV/c
	return momentumHadron;
}

double PhotoPionProduction::crossection(double eps, bool onProton) const {
	const double m = mass(onProton);
	const double s = m * m + 2. * m * eps;
	if (s < sMin())
		return 0.;
	double cross_res = 0.;
	double cross_dir = 0.;
	double cross_dir1 = 0.;
	double cross_dir2 = 0.;
	double sig_res[9];

	// first half of array: 9x proton resonance data | second half of array 9x neutron resonance data
	static const double AMRES[18] = {1.231, 1.440, 1.515, 1.525, 1.675, 1.680, 1.690, 1.895, 1.950, 1.231, 1.440, 1.515, 1.525, 1.675, 1.675, 1.690, 1.895, 1.950};
	static const double BGAMMA[18] = {5.6, 0.5, 4.6, 2.5, 1.0, 2.1, 2.0, 0.2, 1.0, 6.1, 0.3, 4.0, 2.5, 0.0, 0.2, 2.0, 0.2, 1.0};
	static const double WIDTH[18] = {0.11, 0.35, 0.11, 0.1, 0.16, 0.125, 0.29, 0.35, 0.3, 0.11, 0.35, 0.11, 0.1, 0.16, 0.150, 0.29, 0.35, 0.3};
	static const double RATIOJ[18] = {1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2., 1., 0.5, 1., 0.5, 0.5, 1.5, 1., 1.5, 2.};
	static const double AM2[2] = {0.882792, 0.880351};

	const int idx = onProton? 0 : 9;
	double SIG0[9];
	for (int i = 0; i < 9; ++i) {
		SIG0[i] = 4.893089117 / AM2[int(onProton)] * RATIOJ[i + idx] * BGAMMA[i + idx];
	}
	if (eps <= 10.) {
		cross_res = breitwigner(SIG0[0], WIDTH[0 + idx], AMRES[0 + idx], eps, onProton) * Ef(eps, 0.152, 0.17);
		sig_res[0] = cross_res;
		for (int i = 1; i < 9; ++i) {
			sig_res[i] = breitwigner(SIG0[i], WIDTH[i + idx], AMRES[i + idx], eps, onProton) * Ef(eps, 0.15, 0.38);
			cross_res += sig_res[i];
		}
		// direct channel
		if ((eps > 0.1) && (eps < 0.6)) {
			cross_dir1 = 92.7 * Pl(eps, 0.152, 0.25, 2.0)  // single pion production
					   + 40. * std::exp(-(eps - 0.29) * (eps - 0.29) / 0.002)
					   - 15. * std::exp(-(eps - 0.37) * (eps - 0.37) / 0.002);
		} else {
			cross_dir1 = 92.7 * Pl(eps, 0.152, 0.25, 2.0);  // single pion production
		}
		cross_dir2 = 37.7 * Pl(eps, 0.4, 0.6, 2.0);  // double pion production
		cross_dir = cross_dir1 + cross_dir2;
	}
	// fragmentation 2:
	double cross_frag2 = onProton? 80.3 : 60.2;
	cross_frag2 *= Ef(eps, 0.5, 0.1) * std::pow(s, -0.34);
	// multipion production/fragmentation 1 cross section
	double cs_multidiff = 0.;
	double cs_multi = 0.;
	double cross_diffr1 = 0.;
	double cross_diffr2 = 0.;
	double cross_diffr = 0.;
	if (eps > 0.85) {
		double ss1 = (eps - 0.85) / 0.69;
		double ss2 = onProton? 29.3 : 26.4;
		ss2 *= std::pow(s, -0.34) + 59.3 * std::pow(s, 0.095);
		cs_multidiff = (1. - std::exp(-ss1)) * ss2;
		cs_multi = 0.89 * cs_multidiff;
		// diffractive scattering:
		cross_diffr1 = 0.099 * cs_multidiff;
		cross_diffr2 = 0.011 * cs_multidiff;
		cross_diffr = 0.11 * cs_multidiff;
		// **************************************
		ss1 = std::pow(eps - 0.85, 0.75) / 0.64;
		ss2 = 74.1 * std::pow(eps, -0.44) + 62. * std::pow(s, 0.08);
		double cs_tmp = 0.96 * (1. - std::exp(-ss1)) * ss2;
		cross_diffr1 = 0.14 * cs_tmp;
		cross_diffr2 = 0.013 * cs_tmp;
		double cs_delta = cross_frag2 - (cross_diffr1 + cross_diffr2 - cross_diffr);
		if (cs_delta < 0.) {
			cross_frag2 = 0.;
			cs_multi += cs_delta;
		} else {
			cross_frag2 = cs_delta;
		}
		cross_diffr = cross_diffr1 + cross_diffr2;
		cs_multidiff = cs_multi + cross_diffr;
	// in the original SOPHIA code, here is a switch for the return argument.
	// Here, only one case (compare in SOPHIA: NDIR=3) is needed.
	}
	return cross_res + cross_dir + cs_multidiff + cross_frag2;
}

double PhotoPionProduction::Pl(double eps, double epsTh, double epsMax, double alpha) const {
	if (epsTh > eps)
		return 0.;
	const double a = alpha * epsMax / epsTh;
	const double prod1 = std::pow((eps - epsTh) / (epsMax - epsTh), a - alpha);
	const double prod2 = std::pow(eps / epsMax, -a);
	return prod1 * prod2;
}

double PhotoPionProduction::Ef(double eps, double epsTh, double w) const {
	const double wTh = w + epsTh;
	if (eps <= epsTh) {
		return 0.;
	} else if ((eps > epsTh) && (eps < wTh)) {
		return (eps - epsTh) / w;
	} else if (eps >= wTh) {
		return 1.;
	} else {
		throw std::runtime_error("error in function Ef");
	}
}

double PhotoPionProduction::breitwigner(double sigma0, double gamma, double DMM, double epsPrime, bool onProton) const {
	const double m = mass(onProton);
	const double s = m * m + 2. * m * epsPrime;
	const double gam2s = gamma * gamma * s;
	return sigma0 * (s / epsPrime / epsPrime) * gam2s / ((s - DMM * DMM) * (s - DMM * DMM) + gam2s);
}

double PhotoPionProduction::functs(double s, bool onProton) const {
	const double m = mass(onProton);
	const double factor = s - m * m;
	const double epsPrime = factor / 2. / m;
	const double sigmaPg = crossection(epsPrime, onProton);
	return factor * sigmaPg;
}

double PhotoPionProduction::mass(bool onProton) const {
	const double m =  onProton ? mass_proton : mass_neutron;
	return m / GeV * c_squared;
}

double PhotoPionProduction::sMin() const {
	return 1.1646; // [GeV^2] head-on collision
}

void PhotoPionProduction::setSampleLog(bool b) {
	sampleLog = b;
}

void PhotoPionProduction::setCorrectionFactor(double factor) {
	correctionFactor = factor;
}

ref_ptr<PhotonField> PhotoPionProduction::getPhotonField() const {
	return photonField;
}

bool PhotoPionProduction::getHavePhotons() const {
	return havePhotons;
}

bool PhotoPionProduction::getHaveNeutrinos() const {
	return haveNeutrinos;
}

bool PhotoPionProduction::getHaveElectrons() const {
	return haveElectrons;
}

bool PhotoPionProduction::getHaveAntiNucleons() const {
	return haveAntiNucleons;
}

bool PhotoPionProduction::getHaveRedshiftDependence() const {
	return haveRedshiftDependence;
}

double PhotoPionProduction::getLimit() const {
	return limit;
}

bool PhotoPionProduction::getSampleLog() const {
	return sampleLog;
}

double PhotoPionProduction::getCorrectionFactor() const {
	return correctionFactor;
}

void PhotoPionProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string PhotoPionProduction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa

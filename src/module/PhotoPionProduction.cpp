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

PhotoPionProduction::PhotoPionProduction(PhotonField field, bool photons, bool neutrinos, bool electrons, bool antiNucleons, double l, bool redshift) {
	havePhotons = photons;
	haveNeutrinos = neutrinos;
	haveElectrons = electrons;
	haveAntiNucleons = antiNucleons;
	haveRedshiftDependence = redshift;
	limit = l;
	setPhotonField(field);
}

void PhotoPionProduction::setPhotonField(PhotonField field) {
	photonField = field;
	if (haveRedshiftDependence) {
		if (photonField.getHasRedshiftDependence()){
			std::cout << "PhotoPionProduction: tabulated redshift dependence not needed for CMB, switching off" << std::endl;
			haveRedshiftDependence = false;
		}
		else {
			KISS_LOG_WARNING << "PhotoPionProduction: You are using the 2-dimensional tabulated redshift evolution, which is not available for other interactions. To be consistent across all interactions you may deactivate this <setHaveRedshiftDependence(False)>.";
		}
	}
	std::string fname = photonField.getFieldName();
	setDescription("PhotoPionProduction: " + fname);
	if (haveRedshiftDependence)
		initRate(getDataPath("PhotoPionProduction/rate_" + fname.replace(0, 3, "IRBz") + ".txt"));
	else
		initRate(getDataPath("PhotoPionProduction/rate_" + fname + ".txt"));

	int background = (fname == "CMB") ? 1 : 2; // photon background: 1 for CMB, 2 for Kneiske IRB
	this->photonFieldSampling = PhotonFieldSampling(background);
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
		rate = interpolate(gamma, tabLorentz, tabRate) * photonField.getRedshiftScaling(z);

	// cosmological scaling
	rate *= pow(1 + z, 2);

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
	double E_threshold = (photonField.getFieldName() == "CMB") ? 3.72e18 * eV : 5.83e15 * eV;
	if (EpA * (1 + z) < E_threshold)
		return;

	// SOPHIA - input:
	int nature = 1 - static_cast<int>(onProton);  // 0=proton, 1=neutron
	double Ein = EpA / GeV;  // GeV is the SOPHIA standard unit
	double eps = photonFieldSampling.sample_eps(onProton, Ein, z) / GeV;  // GeV for SOPHIA

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
					candidate->addSecondary(-sign * nucleusId(1, 14 + pType), Eout, pos);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (anti-nucleon production)\n" << "Something went wrong in the PhotoPionProduction\n"<< "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			break;
		case 1: // photon
			if (havePhotons)
				candidate->addSecondary(22, Eout, pos);
			break;
		case 2: // positron
			if (haveElectrons)
				candidate->addSecondary(sign * -11, Eout, pos);
			break;
		case 3: // electron
			if (haveElectrons)
				candidate->addSecondary(sign * 11, Eout, pos);
			break;
		case 15: // nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * 12, Eout, pos);
			break;
		case 16: // anti-nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * -12, Eout, pos);
			break;
		case 17: // nu_mu
			if (haveNeutrinos)
				candidate->addSecondary(sign * 14, Eout, pos);
			break;
		case 18: // anti-nu_mu
			if (haveNeutrinos)
				candidate->addSecondary(sign * -14, Eout, pos);
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
					candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (primary particle, A!=1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			}
		} else {  // nucleon is secondary proton or neutron
			candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos);
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

} // namespace crpropa

#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "crpropa/PhotonBackground.h"

#include <kiss/convert.h>
#include <kiss/logger.h>
#include "sophia.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

namespace crpropa {

PhotoPionProduction::PhotoPionProduction(PhotonField field,
		bool photons, bool neutrinos, bool electrons, bool antiNucleons, double l) {
	havePhotons = photons;
	haveNeutrinos = neutrinos;
	haveElectrons = electrons;
	haveAntiNucleons = antiNucleons;
	limit = l;
	setPhotonField(field);
	this->customPhotonField = CustomPhotonField(getDataPath("Scaling/" + photonFieldName(field) + ".txt"));
}

void PhotoPionProduction::setPhotonField(PhotonField field) {
	photonField = field;
	std::string fname = photonFieldName(field);
	setDescription("PhotoPionProduction: " + fname);
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

	infile.close();
}

double PhotoPionProduction::nucleonMFP(double gamma, double z, bool onProton) const {
	const std::vector<double> &tabRate = (onProton)? tabProtonRate : tabNeutronRate;

	// scale nucleus energy instead of background photon energy
	gamma *= (1 + z);
	if (gamma < tabLorentz.front() or (gamma > tabLorentz.back()))
		return std::numeric_limits<double>::max();

	double rate = interpolate2d(z, gamma, tabRedshifts, tabLorentz, tabRate);

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
		if (meanFreePath == std::numeric_limits<double>::max())
			return;
		
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

	// SOPHIA simulates interactions only for protons / neutrons
	// for anti-protons / neutrons assume charge symmetry and change all
	// interaction products from particle <--> anti-particle
	int sign = (id > 0) ? 1 : -1;

	// SOPHIA - input:
	int nature = 1 - static_cast<int>(onProton);  // 0=proton, 1=neutron
	double Ein = EpA / GeV;
	double eps = customPhotonField.sampleEps(onProton, EpA, z);

	// SOPHIA - output:
	double outputEnergy[2000];
	int outPartID[2000];
	int nOutPart;

	#pragma omp critical
	{
		sophiaevent_(nature, Ein, eps, outputEnergy, outPartID, nOutPart);
	}

	// output particle treatment
	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	bool isPrimary = true;  // the first hadron in output is declared the primary
	for (int i = 0; i < nOutPart; i++) { // loop over out-going particles
		double Eout = outputEnergy[i] * GeV;  // only the energy is used; could be changed for more detail
		int pType = outPartID[i];
		try
		{
			switch (pType) {
				case 13: // proton
				case 14: // neutron
					if (isPrimary) {
						if (A == 1) {
							// single interacting nucleon
							candidate->current.setEnergy(Eout);
							candidate->current.setId(sign * nucleusId(1, 14 - pType));
						} else {
							// interacting nucleon is part of nucleus: it is emitted from the nucleus
							candidate->current.setEnergy(E - EpA);
							candidate->current.setId(sign * nucleusId(A - 1, Z - int(onProton)));
							candidate->addSecondary(sign * nucleusId(1, 14 - pType), Eout, pos);
						}
						isPrimary = false;
					} else { 
						candidate->addSecondary(sign * nucleusId(1, 14 - pType), Eout, pos);
					}
					break;
				case -13: // anti-proton
					if (haveAntiNucleons)
						candidate->addSecondary(-sign * nucleusId(1, 1), Eout, pos);
					break;
				case -14: // anti-neutron
					if (haveAntiNucleons)
						candidate->addSecondary(-sign * nucleusId(1, 0), Eout, pos);
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
				case 16: // antinu_e
					if (haveNeutrinos)
						candidate->addSecondary(sign * -12, Eout, pos);
					break;
				case 17: // nu_muon
					if (haveNeutrinos)
						candidate->addSecondary(sign * 14, Eout, pos);
					break;
				case 18: // antinu_muon
					if (haveNeutrinos)
						candidate->addSecondary(sign * -14, Eout, pos);
					break;
				default:
					throw std::runtime_error("PhotoPionProduction: unexpected particle " + kiss::str(pType));
			}
		}
		catch (std::runtime_error &e)
		{
			KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
			throw;
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

} // namespace crpropa

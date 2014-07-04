#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"

#include <kiss/convert.h>
#include "sophia.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {

PhotoPionProduction::PhotoPionProduction(PhotonField field, bool photons,
		bool neutrinos, bool antiNucleons, double l) {
	photonField = field;
	havePhotons = photons;
	haveNeutrinos = neutrinos;
	haveAntiNucleons = antiNucleons;
	limit = l;
	init();
}

void PhotoPionProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	init();
}

void PhotoPionProduction::setHavePhotons(bool b) {
	havePhotons = b;
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

void PhotoPionProduction::init() {
	if (photonField == CMB) {
		init(getDataPath("photopion_CMB.txt"));
		setDescription("PhotoPionProduction: CMB");
	} else if (photonField == IRB) {
		init(getDataPath("photopion_KneiskeIRB.txt"));
		setDescription("PhotoPionProduction: IRB");
	} else {
		throw std::runtime_error(
				"PhotoPionProduction: unknown photon background");
	}
}

void PhotoPionProduction::init(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"PhotoPionProduction: could not open file " + filename);

	// clear previously loaded tables
	tabLorentz.clear();
	tabProtonRate.clear();
	tabNeutronRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b, c;
			infile >> a >> b >> c;
			if (infile) {
				tabLorentz.push_back(pow(10, a));
				tabProtonRate.push_back(b / Mpc);
				tabNeutronRate.push_back(c / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

double PhotoPionProduction::nucleiModification(int A, int X) const {
	if (A == 1)
		return 1.;
	if (A <= 8)
		return 0.85 * pow(X, 2. / 3.);
	return 0.85 * X;
}

void PhotoPionProduction::process(Candidate *candidate) const {
	// the loop should be processed at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (not (isNucleus(id)))
			return;

		double z = candidate->getRedshift();
		double gamma = (1 + z) * candidate->current.getLorentzFactor();

		// check if in tabulated energy range
		if (gamma < tabLorentz.front() or (gamma > tabLorentz.back()))
			return;

		// find interaction with minimum random distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();
		int channel; // interacting particle: 1 for proton, 0 for neutron
		double totalRate = 0;

		// comological scaling of interaction distance (comoving)
		double scaling = pow(1 + z, 2) * photonFieldScaling(photonField, z);

		int A = massNumber(id);
		int Z = chargeNumber(id);
		int N = A - Z;

		// check for interaction on protons
		if (Z > 0) {
			double rate = interpolate(gamma, tabLorentz, tabProtonRate);
			rate *= nucleiModification(A, Z);
			rate *= scaling;
			totalRate += rate;
			channel = 1;
			randDistance = -log(random.rand()) / rate;
		}

		// check for interaction on neutrons
		if (N > 0) {
			double rate = interpolate(gamma, tabLorentz, tabNeutronRate);
			rate *= nucleiModification(A, N);
			rate *= scaling;
			totalRate += rate;
			double d = -log(random.rand()) / rate;
			if (d < randDistance) {
				randDistance = d;
				channel = 0;
			}
		}

		// check if interaction does not happen
		if (step < randDistance) {
			candidate->limitNextStep(limit / totalRate);
			return;
		}

		// interact and repeat with remaining step
		performInteraction(candidate, channel);
		step -= randDistance;
	} while (step > 0);
}

void PhotoPionProduction::performInteraction(Candidate *candidate,
		int channel) const {

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

	// arguments for sophia
	int nature = 1 - channel; // interacting particle: 0 for proton, 1 for neutron
	double Ein = EpA / GeV; // energy of in-going nucleon in GeV
	double momentaList[5][2000]; // momentum list, what are the five components?
	int particleList[2000]; // particle id list
	int nParticles; // number of outgoing particles
	double maxRedshift = 100; // IR photon density is zero above this redshift
	int dummy1; // not needed
	double dummy2[2]; // not needed
	int background = (photonField == CMB) ? 1 : 2; // photon background: 1 for CMB, 2 for Kneiske IRB

#pragma omp critical
	{
		sophiaevent_(nature, Ein, momentaList, particleList, nParticles, z,
				background, maxRedshift, dummy1, dummy2, dummy2);
	}

	for (int i = 0; i < nParticles; i++) { // loop over out-going particles
		double Eout = momentaList[3][i] * GeV; // only the energy is used; could be changed for more detail
		int pType = particleList[i];
		switch (pType) {
		case 13: // proton
		case 14: // neutron
			if (A == 1) {
				// single interacting nucleon
				candidate->current.setEnergy(Eout);
				candidate->current.setId(sign * nucleusId(1, 14 - pType));
			} else {
				// interacting nucleon is part of nucleus: it is emitted from the nucleus
				candidate->current.setEnergy(E - EpA);
				candidate->current.setId(sign * nucleusId(A - 1, Z - channel));
				candidate->addSecondary(sign * nucleusId(1, 14 - pType), Eout);
			}
			break;
		case -13: // anti-proton
		case -14: // anti-neutron
			if (haveAntiNucleons)
				candidate->addSecondary(-sign * nucleusId(1, 14 + pType), Eout);
			break;
		case 1: // photon
			if (havePhotons)
				candidate->addSecondary(22, Eout);
			break;
		case 2: // positron
			if (havePhotons)
				candidate->addSecondary(sign * -11, Eout);
			break;
		case 3: // electron
			if (havePhotons)
				candidate->addSecondary(sign * 11, Eout);
			break;
		case 15: // nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * 12, Eout);
			break;
		case 16: // antinu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * -12, Eout);
			break;
		case 17: // nu_muon
			if (haveNeutrinos)
				candidate->addSecondary(sign * 14, Eout);
			break;
		case 18: // antinu_muon
			if (haveNeutrinos)
				candidate->addSecondary(sign * -14, Eout);
			break;
		default:
			throw std::runtime_error(
					"PhotoPionProduction: unexpected particle "
							+ kiss::str(pType));
		}
	}
}

double PhotoPionProduction::lossLength(int id, double gamma, double z) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;

	gamma *= (1 + z); // cosmological scaling of photon energy
	if (gamma < tabLorentz.front() or (gamma > tabLorentz.back()))
		return std::numeric_limits<double>::max();

	double lossRate = 0;
	if (Z > 0)
		lossRate += interpolate(gamma, tabLorentz, tabProtonRate)
				* nucleiModification(A, Z);
	if (N > 0)
		lossRate += interpolate(gamma, tabLorentz, tabProtonRate)
				* nucleiModification(A, N);

	// protons / neutrons keep as energy the fraction of mass to delta-resonance mass
	// nuclei approximately lose the energy that the interacting nucleon is carrying
	double relativeEnergyLoss = (A == 1) ? 1 - 938. / 1232. : 1. / A;
	lossRate *= relativeEnergyLoss;

	// cosmological scaling of photon density
	lossRate *= pow(1 + z, 3) * photonFieldScaling(photonField, z);

	return 1. / lossRate;
}

} // namespace crpropa

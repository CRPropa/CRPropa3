#include "mpc/module/PhotoPionProduction.h"

#include <kiss/convert.h>
#include "sophia.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace mpc {

PhotoPionProduction::PhotoPionProduction(PhotonField photonField) :
		photonField(photonField) {
	init();
}

void PhotoPionProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	init();
}

void PhotoPionProduction::init() {
	if (photonField == CMB) {
		init(getDataPath("photopion_CMB.txt"));
		setDescription("PhotoPionProduction: CMB");
	}
	else if (photonField == IRB) {
		init(getDataPath("photopion_IRB.txt"));
		setDescription("PhotoPionProduction: IRB");
	}
	else
		throw std::runtime_error(
				"PhotoPionProduction: unknown photon background");
}

void PhotoPionProduction::init(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"PhotoPionProduction: could not open file " + filename);

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b, c;
			infile >> a >> b >> c;
			if (infile) {
				energy.push_back(a * EeV);
				pRate.push_back(b / Mpc);
				nRate.push_back(c / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	infile.close();
}

bool PhotoPionProduction::setNextInteraction(Candidate *candidate,
		InteractionState &interaction) const {
	if (not(candidate->current.isNucleus()))
		return false; // accept only nuclei

	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	int N = A - Z;
	double EpA = E / A * (1 + z); // CMB energies increase with (1+z)^3

	// check if out of energy range
	if ((EpA < energy.front()) or (EpA > energy.back()))
		return false;

	// find interaction with minimum random distance
	interaction.distance = std::numeric_limits<double>::max();
	Random &random = Random::instance();

	// check for interaction on protons
	if (Z > 0) {
		double rate = interpolate(EpA, energy, pRate);
		if (rate > 0) {
			if (A > 1) {
				if (A < 8)
					rate *= 0.85 * pow(Z, 2. / 3.);
				if (A >= 8)
					rate *= 0.85 * Z;
			}
			interaction.distance = -log(random.rand()) / rate;
			interaction.channel = 1;
		}
	}

	// check for interaction on neutrons
	if (N > 0) {
		double rate = interpolate(EpA, energy, nRate);
		if (rate > 0) {
			if (A > 1) {
				if (A < 8)
					rate *= 0.85 * pow(N, 2. / 3.);
				if (A >= 8)
					rate *= 0.85 * N;
			}

			double d = -log(random.rand()) / rate;
			if (d < interaction.distance) {
				interaction.distance = -log(random.rand()) / rate;
				interaction.channel = 0;
			}
		}
	}

	// interaction length is proportional to 1 / (photon density)
	interaction.distance /= photonFieldScaling(photonField, z);
	// convert to comoving frame
	interaction.distance *= (1 + z);

	candidate->setInteractionState(getDescription(), interaction);
	return true;
}

void PhotoPionProduction::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(getDescription(), interaction);
	candidate->clearInteractionStates();

	// charge number loss of interaction nucleus
	int dZ = interaction.channel;
	// final proton number of emitted nucleon
	int Zfinal = dZ;
	// 50% probability of isospin change p <-> n
	Random &random = Random::instance();
	if (random.rand() < 1. / 2.)
		Zfinal = abs(Zfinal - 1);

	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();

	if (A == 1) {
		// interaction on single nucleon
		candidate->current.setEnergy(E * 938. / 1232.);
		candidate->current.setId(nucleusId(1, Zfinal));
	} else {
		// interaction on nucleus, update nucleus and emit nucleon
		candidate->current.setEnergy(E * (A - 1) / A);
		candidate->current.setId(nucleusId(A - 1, Z - dZ));
		candidate->addSecondary(nucleusId(1, Zfinal), E / A * 938. / 1232.);
	}
}

double PhotoPionProduction::energyLossLength(int id, double E) {
	int A = massNumberFromNucleusId(id);
	int Z = chargeNumberFromNucleusId(id);
	int N = A - Z;

	double EpA = E / A;
	if ((EpA < energy.front()) or (EpA > energy.back()))
		return std::numeric_limits<double>::max();

	double lossRate = 0;
	double relativeEnergyLoss = 1. / double(A);

	if (Z > 0) {
		double rate = interpolate(EpA, energy, pRate);
		if (A > 1)
			if (A < 8)
				rate *= 0.85 * pow(Z, 2. / 3.);
			if (A >= 8)
				rate *= 0.85 * Z;
		lossRate += relativeEnergyLoss * rate;
	}

	if (N > 0) {
		double rate = interpolate(EpA, energy, nRate);
		if (A > 1)
			if (A < 8)
				rate *= 0.85 * pow(N, 2. / 3.);
			if (A >= 8)
				rate *= 0.85 * N;
		lossRate += relativeEnergyLoss * rate;
	}

	return 1. / lossRate;
}

SophiaPhotoPionProduction::SophiaPhotoPionProduction(PhotonField photonField,
		bool photons, bool neutrinos, bool antiNucleons) :
		PhotoPionProduction(photonField), havePhotons(photons), haveNeutrinos(
				neutrinos), haveAntiNucleons(antiNucleons) {
}

void SophiaPhotoPionProduction::setHavePhotons(bool b) {
	havePhotons = b;
}

void SophiaPhotoPionProduction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void SophiaPhotoPionProduction::setHaveAntiNucleons(bool b) {
	haveAntiNucleons = b;
}

void SophiaPhotoPionProduction::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(getDescription(), interaction);
	candidate->clearInteractionStates();

	int channel = interaction.channel; // 1 for interaction proton, 0 for neutron

	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	double EpA = E / A;

	// arguments for sophia
	int nature = 1 - channel; // interacting particle: 0 for proton, 1 for neutron
	double Ein = EpA / GeV; // energy of in-going nucleon in GeV
	double momentaList[5][2000]; // momentum list, what are the five components?
	int particleList[2000]; // particle id list
	int nParticles; // number of outgoing particles
	double redshift = candidate->getRedshift();

	int background; // Photon background: 1 for CMB, 2 for Kneiske IRB
	if (photonField == CMB)
		background = 1;
	else if (photonField == IRB)
		background = 2;
	else
		throw std::runtime_error(
				"SophiaPhotoPionProduction: Only CMB an IRB provided");

	double maxRedshift = 100; // IR photon density is zero above this redshift
	int dummy1;
	double dummy2[2];

#pragma omp critical
	{
		sophiaevent_(nature, Ein, momentaList, particleList, nParticles,
				redshift, background, maxRedshift, dummy1, dummy2, dummy2);
	}

	for (int i = 0; i < nParticles; i++) { // loop over out-going particles
		double Eout = momentaList[3][i] * GeV; // only the energy is used; could be changed for more detail
		int pType = particleList[i];

		switch (pType) {
		case 13: // proton
		case 14: // neutron
			if (A == 1) { // in-going particle was a nucleon: update its properties
				candidate->current.setEnergy(Eout);
				candidate->current.setId(nucleusId(1, 14 - pType));
			} else { // in-going particle was a nucleus: update nucleus and emit nucleon
				candidate->current.setEnergy(E - Eout);
				candidate->current.setId(nucleusId(A - 1, Z - channel));
				candidate->addSecondary(nucleusId(1, 14 - pType), Eout);
			}
			break;
		case -13: // anti-proton
			if (haveAntiNucleons)
				candidate->addSecondary(-nucleusId(1, 1), Eout);
			break;
		case -14: // anti-neutron
			if (haveAntiNucleons)
				candidate->addSecondary(-nucleusId(1, 0), Eout);
			break;
		case 1: // photon
			if (havePhotons)
				candidate->addSecondary(22, Eout);
			break;
		case 2: // positron
			if (havePhotons)
				candidate->addSecondary(-11, Eout);
			break;
		case 3: // electron
			if (havePhotons)
				candidate->addSecondary(11, Eout);
			break;
		case 15: // nu_e
			if (haveNeutrinos)
				candidate->addSecondary(12, Eout);
			break;
		case 16: // antinu_e
			if (haveNeutrinos)
				candidate->addSecondary(-12, Eout);
			break;
		case 17: // nu_muon
			if (haveNeutrinos)
				candidate->addSecondary(14, Eout);
			break;
		case 18: // antinu_muon
			if (haveNeutrinos)
				candidate->addSecondary(-14, Eout);
			break;
		default:
			throw std::runtime_error(
					"PhotoPionProduction: unexpected particle "
							+ kiss::str(pType));
		}
	}
}

} // namespace mpc

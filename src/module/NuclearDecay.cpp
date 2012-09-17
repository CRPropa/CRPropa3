#include "mpc/module/NuclearDecay.h"
#include "mpc/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace mpc {

NuclearDecay::NuclearDecay(bool electrons, bool neutrinos) :
		haveElectrons(electrons), haveNeutrinos(neutrinos) {
	setDescription("NuclearDecay");

	// load decay table
	std::string filename = getDataPath("/NuclearDecay/decayTable.txt");
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"mpc::NuclearDecay: could not open file " + filename);

	decayTable.resize(27 * 31);
	while (infile.good()) {
		if (infile.peek() != '#') {
			InteractionState decay;
			int Z, N;
			infile >> Z >> N >> decay.channel >> decay.distance;
			decay.distance *= c_light; // mean decay distance [m]
			if (infile)
				decayTable[Z * 31 + N].push_back(decay);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();

	// generate inverse cdf for electron kinetic energies in neutron decays
	double Q = mass_neutron - mass_proton;
	double cdf = 0;
	for (int i = 0; i < 50; i++) {
		double E = mass_electron + i / 50. * (Q - mass_electron);
		cdf += sqrt(pow(E, 2) - pow(mass_electron, 2)) * pow(Q - E, 2) * E;
		cdfBeta.push_back(cdf);
		tBeta.push_back((E - mass_electron) * c_squared);
	}
	for (int i = 0; i < 50; i++) {
		cdfBeta[i] /= cdf;
	}
}

void NuclearDecay::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void NuclearDecay::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

bool NuclearDecay::setNextInteraction(Candidate *candidate,
		InteractionState &interaction) const {
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	int N = A - Z;

	std::vector<InteractionState> decays = decayTable[Z * 31 + N];
	if (decays.size() == 0)
		return false;

	// find interaction mode with minimum random decay distance
	Random &random = Random::instance();
	interaction.distance = std::numeric_limits<double>::max();
	for (size_t i = 0; i < decays.size(); i++) {
		double d = -log(random.rand()) * decays[i].distance;
		if (d > interaction.distance)
			continue;
		interaction.distance = d;
		interaction.channel = decays[i].channel;
	}

	// special relativistic time dilation
	interaction.distance *= candidate->current.getLorentzFactor();
	// convert to comoving frame
	interaction.distance *= (1 + candidate->getRedshift());

	candidate->setInteractionState(getDescription(), interaction);
	return true;
}

void NuclearDecay::performInteraction(Candidate *candidate) const {
	InteractionState decay;
	candidate->getInteractionState(getDescription(), decay);
	candidate->clearInteractionStates();

	// parse decay channels
	int nBetaMinus = digit(decay.channel, 10000);
	int nBetaPlus = digit(decay.channel, 1000);
	int nAlpha = digit(decay.channel, 100);
	int nProton = digit(decay.channel, 10);
	int nNeutron = digit(decay.channel, 1);

	// perform decays
	for (size_t i = 0; i < nBetaMinus; i++)
		betaDecay(candidate, false);
	for (size_t i = 0; i < nBetaPlus; i++)
		betaDecay(candidate, true);
	for (size_t i = 0; i < nAlpha; i++)
		nucleonEmission(candidate, 4, 2);
	for (size_t i = 0; i < nProton; i++)
		nucleonEmission(candidate, 1, 1);
	for (size_t i = 0; i < nNeutron; i++)
		nucleonEmission(candidate, 1, 0);
}

void NuclearDecay::betaDecay(Candidate *candidate, bool isBetaPlus) const {
	double gamma = candidate->current.getLorentzFactor();
	int Z = candidate->current.getChargeNumber();
	int A = candidate->current.getMassNumber();
	double mass = candidate->current.getMass();

	// beta- decay
	int electronId = 11; // electron
	int neutrinoId = -12; // anti-electron neutrino
	int dZ = 1;
	// beta+ decay
	if (isBetaPlus) {
		electronId = -11; // positron
		neutrinoId = 12; // electron neutrion
		dZ = -1;
	}

	// update candidate
	candidate->current.setId(getNucleusId(A, Z + dZ));
	candidate->current.setLorentzFactor(gamma);

	// random kinetic energy of electron in neutron decay
	Random &random = Random::instance();
	double T = interpolate(random.rand(), cdfBeta, tBeta);
	double Q = (mass - candidate->current.getMass() - mass_electron)
			* c_squared;
	double Qneutron = (mass_neutron - mass_proton - mass_electron) * c_squared;
	// electron energy in this decay
	double E = T * Q / Qneutron + mass_electron * c_squared;
	double p = sqrt(E * E - pow(mass_electron * c_squared, 2));
	double cosTheta = 2 * random.rand() - 1;

	if (haveElectrons)
		// add electron/positron boosted to lab frame
		candidate->addSecondary(electronId, gamma * E * (1 + p * cosTheta));
	if (haveNeutrinos)
		// add neutrino with remaining energy and opposite momentum
		candidate->addSecondary(neutrinoId,
				gamma * (Q - E) * (1 - p * cosTheta));
}

void NuclearDecay::nucleonEmission(Candidate *candidate, int dA, int dZ) const {
	int Z = candidate->current.getChargeNumber();
	int A = candidate->current.getMassNumber();
	double EpA = candidate->current.getEnergy()
			/ double(candidate->current.getMassNumber());
	candidate->current.setId(getNucleusId(A - dA, Z - dZ));
	candidate->current.setEnergy(EpA * (A - dA));
	candidate->addSecondary(getNucleusId(dA, dZ), EpA * dA);
}

} // namespace mpc

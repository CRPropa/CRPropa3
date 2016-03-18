#include "crpropa/module/NuclearDecay.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {

NuclearDecay::NuclearDecay(bool electrons, bool neutrinos, double l) {
	haveElectrons = electrons;
	haveNeutrinos = neutrinos;
	limit = l;
	setDescription("NuclearDecay");

	// load decay table
	std::string filename = getDataPath("nuclear_decay.txt");
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"crpropa::NuclearDecay: could not open file " + filename);

	decayTable.resize(27 * 31);
	while (infile.good()) {
		if (infile.peek() != '#') {
			DecayMode decay;
			int Z, N;
			double lifetime;
			infile >> Z >> N >> decay.channel >> lifetime;
			decay.rate = 1. / lifetime / c_light; // decay rate in [1/m]
			if (infile)
				decayTable[Z * 31 + N].push_back(decay);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

void NuclearDecay::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void NuclearDecay::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void NuclearDecay::setLimit(double l) {
	limit = l;
}

void NuclearDecay::process(Candidate *candidate) const {
	// the loop should be processed at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (not (isNucleus(id)))
			return;

		int A = massNumber(id);
		int Z = chargeNumber(id);
		int N = A - Z;

		// check if particle can decay
		std::vector<DecayMode> decays = decayTable[Z * 31 + N];
		if (decays.size() == 0)
			return;

		// find interaction mode with minimum random decay distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();
		int channel;
		double totalRate = 0;

		for (size_t i = 0; i < decays.size(); i++) {
			double rate = decays[i].rate;
			rate /= candidate->current.getLorentzFactor();
			totalRate += rate;
			double d = -log(random.rand()) / rate;
			if (d > randDistance)
				continue;
			randDistance = d;
			channel = decays[i].channel;
		}

		// check if interaction doesn't happen
		if (step < randDistance) {
			// limit next step to a fraction of the mean free path
			candidate->limitNextStep(limit / totalRate);
			return;
		}

		// interact and repeat with remaining step
		performInteraction(candidate, channel);
		step -= randDistance;
	} while (step > 0);
}

void NuclearDecay::performInteraction(Candidate *candidate, int channel) const {
	// interpret decay channel
	int nBetaMinus = digit(channel, 10000);
	int nBetaPlus = digit(channel, 1000);
	int nAlpha = digit(channel, 100);
	int nProton = digit(channel, 10);
	int nNeutron = digit(channel, 1);

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
	int id = candidate->current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);
	double mass = candidate->current.getMass();

	// beta- decay
	int electronId = 11; // electron
	int neutrinoId = -12; // anti-electron neutrino
	int dZ = 1;
	// beta+ decay
	if (isBetaPlus) {
		electronId = -11; // positron
		neutrinoId = 12; // electron neutrino
		dZ = -1;
	}

	// update candidate, nuclear recoil negligible
	candidate->current.setId(nucleusId(A, Z + dZ));
	candidate->current.setLorentzFactor(gamma);

	if (not (haveElectrons or haveNeutrinos))
		return;

	// Q-value of the decay
	double newMass = candidate->current.getMass();
	double Q = (mass - newMass - mass_electron) * c_squared;

	// generate cdf of electron energy, neglecting Coulomb correction
	// see Basdevant, Fundamentals in Nuclear Physics, eq. (4.92)
	std::vector<double> energies;
	std::vector<double> densities; // cdf(E), unnormalized

	double me = mass_electron * c_squared;
	double cdf = 0;
	for (int i = 0; i <= 50; i++) {
		double E = me + i / 50. * Q;
		cdf += E * sqrt(E * E - me * me) * pow(Q + me - E, 2);
		energies.push_back(E);
		densities.push_back(cdf);
	}

	// draw random electron energy and angle
	Random &random = Random::instance();
	double E = interpolate(random.rand() * cdf, densities, energies);
	double p = sqrt(E * E - me * me);  // p*c
	double cosTheta = 2 * random.rand() - 1;

	// boost to lab frame
	double Ee = gamma * (E - p * cosTheta);
	double Enu = gamma * (Q + me - E) * (1 + cosTheta);  // pnu*c ~ Enu

  Vector3d pos = randomPositionInPropagationStep(candidate);
	if (haveElectrons)
		candidate->addSecondary(electronId, Ee, pos);
	if (haveNeutrinos)
		candidate->addSecondary(neutrinoId, Enu, pos);
}

void NuclearDecay::nucleonEmission(Candidate *candidate, int dA, int dZ) const {
	int id = candidate->current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);
	double EpA = candidate->current.getEnergy() / double(A);
	candidate->current.setId(nucleusId(A - dA, Z - dZ));
	candidate->current.setEnergy(EpA * (A - dA));
  Vector3d pos = randomPositionInPropagationStep(candidate);
	candidate->addSecondary(nucleusId(dA, dZ), EpA * dA, pos);
}

} // namespace crpropa

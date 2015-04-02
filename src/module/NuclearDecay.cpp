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

	// generate cdf for electron kinetic energies in free neutron decays
	// see Basdevant, Fundamentals in Nuclear Physics, (4.92)
	double dm = mass_neutron - mass_proton;
	double cdf = 0;
	for (int i = 0; i < 50; i++) {
		double E = mass_electron + i / 50. * (dm - mass_electron);
		cdf += sqrt(pow(E, 2) - pow(mass_electron, 2)) * pow(dm - E, 2) * E;
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
			rate *= candidate->current.getLorentzFactor();
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

	// update candidate, energy loss negligible
	candidate->current.setId(nucleusId(A, Z + dZ));
	candidate->current.setLorentzFactor(gamma);

	if (not (haveElectrons or haveNeutrinos))
		return;

	// random electron kinetic energy in free neutron decay
	Random &random = Random::instance();
	double Tn = interpolate(random.rand(), cdfBeta, tBeta);
	double Qn = (mass_neutron - mass_proton - mass_electron) * c_squared;

	// scale to Q-value of the given nuclear decay (rough approximation)
	double Q = (mass - candidate->current.getMass() - mass_electron)
			* c_squared;
	double T = Tn * Q / Qn;

	// electron energy and momentum in the nuclear rest frame
	double E = T + mass_electron * c_squared;
	double p = sqrt(E * E - pow(mass_electron * c_squared, 2));

	// random angle with respect to the cosmic ray direction
	double cosTheta = 2 * random.rand() - 1;

	if (haveElectrons)
		// add electron/positron boosted to lab frame
		candidate->addSecondary(electronId, gamma * (E + p * cosTheta));
	if (haveNeutrinos)
		// add neutrino with remaining energy and opposite momentum
		candidate->addSecondary(neutrinoId, gamma * ((Q - T) - p * cosTheta));
}

void NuclearDecay::nucleonEmission(Candidate *candidate, int dA, int dZ) const {
	int id = candidate->current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);
	double EpA = candidate->current.getEnergy() / double(A);
	candidate->current.setId(nucleusId(A - dA, Z - dZ));
	candidate->current.setEnergy(EpA * (A - dA));
	candidate->addSecondary(nucleusId(dA, dZ), EpA * dA);
}

} // namespace crpropa

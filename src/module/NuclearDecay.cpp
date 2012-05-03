#include "mpc/module/NuclearDecay.h"

#include <fstream>
#include <limits>
#include <math.h>
#include <stdexcept>

namespace mpc {

NuclearDecay::NuclearDecay(bool electrons, bool neutrinos) {
	setDescription("NuclearDecay");
	haveElectrons = electrons;
	haveNeutrinos = neutrinos;

	// load decay table
	std::string filename = getDataPath("/NuclearDecay/decay_table.txt");
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"mpc::NuclearDecay: could not open file " + filename);

	while (infile.good()) {
		if (infile.peek() != '#') {
			InteractionState decay;
			int Z, N;
			infile >> Z >> N >> decay.distance >> decay.channel;
			if (infile) {
				decay.distance *= c_light;
				decayTable[getNucleusId(Z + N, Z)].push_back(decay);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();

	// generate inverse cdf for electron kinetic energies in neutron decays
	double Q = mass_neutron - mass_proton;
	double cdf = 0;
	std::vector<double> x, y;
	x.resize(50);
	y.resize(50);
	for (int i = 0; i < 50; i++) {
		double E = mass_electron + i / 50. * (Q - mass_electron);
		cdf += sqrt(pow(E, 2) - pow(mass_electron, 2)) * pow(Q - E, 2) * E;
		x[i] = cdf;
		y[i] = (E - mass_electron) * c_squared;
	}
	for (int i = 0; i < 50; i++) {
		x[i] /= x.back();
	}

	acc = gsl_interp_accel_alloc();
	Tbeta = gsl_spline_alloc(gsl_interp_linear, x.size());
	gsl_spline_init(Tbeta, &x[0], &y[0], x.size());
}

bool NuclearDecay::setNextInteraction(Candidate *candidate,
		InteractionState &decay) const {
	int id = candidate->current.getId();

	std::map<int, std::vector<InteractionState> >::const_iterator iMode =
			decayTable.find(id);
	if (iMode == decayTable.end())
		return false;

	const std::vector<InteractionState> &states = iMode->second;

	// find decay mode with minimum random decay distance
	decay.distance = std::numeric_limits<double>::max();
	int decayChannel;
	for (size_t i = 0; i < states.size(); i++) {
		double d = -log(Random::instance().rand()) * states[i].distance;
		if (d > decay.distance)
			continue;
		decay.distance = d;
		decay.channel = states[i].channel;
	}
	decay.distance *= candidate->current.getLorentzFactor();
	candidate->setInteractionState(getDescription(), decay);
	return true;
}

void NuclearDecay::performInteraction(Candidate *candidate) const {
	InteractionState decay;
	candidate->getInteractionState(getDescription(), decay);
	candidate->clearInteractionStates();

	// parse decay channels
	int nBeta = digit(decay.channel, 10000);
	int nBetaPlus = digit(decay.channel, 1000);
	int nAlpha = digit(decay.channel, 100);
	int nProton = digit(decay.channel, 10);
	int nNeutron = digit(decay.channel, 1);

	// perform decays
	for (size_t i = 0; i < nBeta; i++)
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

	int electronId = 11; // electron
	int neutrinoId = -12; // anti-electron neutrino
	int dZ = 1;
	if (isBetaPlus) {
		electronId = -11; // positron
		neutrinoId = 12; // electron neutrion
		dZ = -1;
	}

	// update candidate
	candidate->current.setId(getNucleusId(A, Z + dZ));
	candidate->current.setLorentzFactor(gamma);

	// random kinetic energy of electron in neutron decay
	double T = gsl_spline_eval(Tbeta, Random::instance().rand(), acc);
	double Q = (mass - candidate->current.getMass() - mass_electron) * c_squared;
	double Qneutron = (mass_neutron - mass_proton - mass_electron) * c_squared;
	// electron energy in this decay
	double E = T * Q / Qneutron + mass_electron * c_squared;
	double p = sqrt(E * E - pow(mass_electron * c_squared, 2));
	double cosTheta = 2 * Random::instance().rand() - 1;

	if (haveElectrons)
		// add electron/positron boosted to lab frame
		candidate->addSecondary(electronId, gamma * E * (1 + p * cosTheta));
	if (haveNeutrinos)
		// add neutrino with remaining energy and opposite momentum
		candidate->addSecondary(neutrinoId, gamma * (Q - E) * (1 - p * cosTheta));
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

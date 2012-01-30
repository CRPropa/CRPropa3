#include "mpc/module/NuclearDecay.h"

#include <fstream>
#include <limits>
#include <math.h>

namespace mpc {

NuclearDecay::NuclearDecay() {
	cached_id = 0;
	// read decay data
	int id;
	DecayMode mode;
	std::ifstream infile("data/NuclearDecay/decay.txt");
	while (infile.good()) {
		if (infile.peek() != '#') {
			infile >> id >> mode.distance >> mode.channel;
			if (infile) {
				mode.distance *= c_light;
				modeMap[id].push_back(mode);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

std::string NuclearDecay::getDescription() const {
	return "Nuclear Decay";
}

void NuclearDecay::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	int id = candidate->current.getId();
	std::vector<DecayMode> modes = modeMap[id];

	// check if stable
	if (modes.size() == 0)
		return;

	// check if a decay is already cached, if not select one
	if (id != cached_id) {
		cached_id = id;
		// find decay mode with minimum random decay distance
		cached_distance = std::numeric_limits<double>::max();
		for (size_t i = 0; i < modes.size(); i++) {
			double d = -log(mtrand.rand()) * modes[i].distance;
			if (d > cached_distance)
				continue;
			cached_distance = d;
			cached_channel = modes[i].channel;
		}
	}

	// check if life-time is over, reduce life-time
	double gamma = candidate->current.getLorentzFactor();
	double step = candidate->getCurrentStep() / gamma;
	if (cached_distance > step) {
		cached_distance -= step;
		candidate->limitNextStep(cached_distance * gamma);
		return;
	}

	// else life time is over, decay
	cached_id = 0;
	decay(candidate, secondaries);
}

void createSecondary(Candidate *candidate,
		std::vector<Candidate *> &secondaries, int id, double energy) {
	ParticleState initial = candidate->current;
	initial.setEnergy(energy);
	initial.setId(id);
	Candidate secondary;
	secondary.current = initial;
	secondary.initial = initial;
	secondary.setNextStep(candidate->getCurrentStep());
	secondaries.push_back(&secondary);
}

void NuclearDecay::decay(Candidate *candidate, std::vector<Candidate *> &secondaries) {
	int nBeta = cached_channel / 10000;
	int nBetaPlus = (cached_channel % 10000) / 1000;
	int nAlpha = (cached_channel % 1000) / 100;
	int nProton = (cached_channel % 100) / 10;
	int nNeutron = cached_channel % 10;

	int dA = -4 * nAlpha - nProton - nNeutron;
	int dZ = -2 * nAlpha - nProton + nBeta - nBetaPlus;

	int id = candidate->current.getId();
	int A = getMassNumberFromNucleusId(id);
	int Z = getChargeNumberFromNucleusId(id);
	double energyPerNucleon = candidate->current.getEnergy() / double(A);

	candidate->current.setId(getNucleusId(A + dA, Z + dZ));
	candidate->current.setEnergy(energyPerNucleon * (A + dA));

	for (size_t i = 0; i < nBeta; i++) {
		// electron + neutrino
	}
	for (size_t i = 0; i < nBetaPlus; i++) {
		// positron + neutrino
	}
	for (size_t i = 0; i < nAlpha; i++) {
		createSecondary(candidate, secondaries, getNucleusId(4, 2),
				energyPerNucleon * 4);
	}
	for (size_t i = 0; i < nProton; i++) {
		createSecondary(candidate, secondaries, getNucleusId(4, 2),
				energyPerNucleon);
	}
	for (size_t i = 0; i < nNeutron; i++) {
		createSecondary(candidate, secondaries, getNucleusId(4, 2),
				energyPerNucleon);
	}
}

} // namespace mpc

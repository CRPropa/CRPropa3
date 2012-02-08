#include "mpc/module/NuclearDecay.h"

#include <fstream>
#include <limits>
#include <math.h>

namespace mpc {

NuclearDecay::NuclearDecay() {
	std::ifstream infile("data/NuclearDecay/decay.txt");
	while (infile.good()) {
		if (infile.peek() != '#') {
			DecayMode mode;
			int Z, N;
			infile >> Z >> N >> mode.distance >> mode.channel;
			if (infile) {
				mode.distance *= c_light;
				modeMap[getNucleusId(Z + N, Z)].push_back(mode);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

std::string NuclearDecay::getDescription() const {
	return "Nuclear Decay";
}

void NuclearDecay::process(Candidate *candidate) {
	double gamma = candidate->current.getLorentzFactor();
	double step = candidate->getCurrentStep() / gamma;

	while (true) {
		int id = candidate->current.getId();

		// set a new interaction if necessary
		if (id != cached_id) {
			if (setNextInteraction(candidate) == false)
				return;
			cached_id = id;
		}
		// if counter not over, reduce and return
		if (cached_distance > step) {
			cached_distance -= step;
			candidate->limitNextStep(cached_distance * gamma);
			return;
		}
		// counter over: interact
		cached_id = 0;
		step -= cached_distance;
		performInteraction(candidate);
	}
}

bool NuclearDecay::setNextInteraction(Candidate *candidate) {
	int id = candidate->current.getId();

	// check if stable
	if (modeMap[id].size() == 0)
		return false;

	// find decay mode with minimum random decay distance
	cached_distance = std::numeric_limits<double>::max();
	for (size_t i = 0; i < modeMap[id].size(); i++) {
		double d = -log(mtrand.rand()) * modeMap[id][i].distance;
		if (d > cached_distance)
			continue;
		cached_distance = d;
		cached_channel = modeMap[id][i].channel;
	}
	return true;
}

void NuclearDecay::performInteraction(Candidate *candidate) {
	int nBeta = cached_channel / 10000;
	int nBetaPlus = (cached_channel % 10000) / 1000;
	int nAlpha = (cached_channel % 1000) / 100;
	int nProton = (cached_channel % 100) / 10;
	int nNeutron = cached_channel % 10;

	int dA = -4 * nAlpha - nProton - nNeutron;
	int dZ = -2 * nAlpha - nProton + nBeta - nBetaPlus;

	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	double EpA = candidate->current.getEnergy() / double(A);

	// update particle
	candidate->current.setId(getNucleusId(A + dA, Z + dZ));
	candidate->current.setEnergy(EpA * (A + dA));

	// create secondaries
	for (size_t i = 0; i < nBeta; i++) {
		// electron + neutrino
	}
	for (size_t i = 0; i < nBetaPlus; i++) {
		// positron + neutrino
	}
	for (size_t i = 0; i < nAlpha; i++) {
		candidate->addSecondary(getNucleusId(4, 2), EpA * 4);
	}
	for (size_t i = 0; i < nProton; i++) {
		candidate->addSecondary(getNucleusId(4, 2), EpA);
	}
	for (size_t i = 0; i < nNeutron; i++) {
		candidate->addSecondary(getNucleusId(4, 2), EpA);
	}
}

} // namespace mpc

#include "mpc/module/NuclearDecay.h"

#include <fstream>
#include <limits>
#include <math.h>

namespace mpc {

NuclearDecay::NuclearDecay() {
	name = "mpc::NuclearDecay";

	std::ifstream infile("data/NuclearDecay/decay.txt");
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
}

std::string NuclearDecay::getDescription() const {
	return name;
}

void NuclearDecay::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	double gamma = candidate->current.getLorentzFactor();
	double step = candidate->getCurrentStep() / gamma;
	InteractionState decay;

	while (true) {
		// set a new decay if necessary
		if (not (candidate->getInteractionState(name, decay))) {
			if (not (setNextInteraction(candidate)))
				return;
		}

		// get the updated decay
		candidate->getInteractionState(name, decay);

		// if not over, reduce distance and return
		if (decay.distance > step) {
			decay.distance -= step;
			candidate->limitNextStep(decay.distance * gamma);
			candidate->setInteractionState(name, decay);
			return;
		}

		// counter over: interact
		step -= decay.distance;
		performInteraction(candidate, secondaries);
	}
}

bool NuclearDecay::setNextInteraction(Candidate *candidate) {
	int id = candidate->current.getId();

	// check if stable
	if (decayTable[id].size() == 0)
		return false;

	// find decay mode with minimum random decay distance
	InteractionState decay;
	decay.distance = std::numeric_limits<double>::max();
	int decayChannel;
	for (size_t i = 0; i < decayTable[id].size(); i++) {
		double d = -log(random.rand()) * decayTable[id][i].distance;
		if (d > decay.distance)
			continue;
		decay.distance = d;
		decay.channel = decayTable[id][i].channel;
	}
	candidate->setInteractionState(name, decay);
	return true;
}

void NuclearDecay::performInteraction(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	InteractionState decay;
	candidate->getInteractionState(name, decay);
	candidate->clearInteractionStates();

	// parse decay channel
	int nBeta = decay.channel / 10000;
	int nBetaPlus = (decay.channel % 10000) / 1000;
	int nAlpha = (decay.channel % 1000) / 100;
	int nProton = (decay.channel % 100) / 10;
	int nNeutron = decay.channel % 10;

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
		// electron + neutrino not implemented
	}
	for (size_t i = 0; i < nBetaPlus; i++) {
		// positron + neutrino not implemented
	}
	for (size_t i = 0; i < nAlpha; i++) {
		addSecondary(secondaries, candidate, getNucleusId(4, 2), EpA * 4);
	}
	for (size_t i = 0; i < nProton; i++) {
		addSecondary(secondaries, candidate, getNucleusId(4, 2), EpA);
	}
	for (size_t i = 0; i < nNeutron; i++) {
		addSecondary(secondaries, candidate, getNucleusId(4, 2), EpA);
	}
}

} // namespace mpc

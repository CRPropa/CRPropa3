#include "mpc/module/NuclearDecay.h"

#include <fstream>
#include <limits>
#include <math.h>
#include <stdexcept>

namespace mpc {

NuclearDecay::NuclearDecay() : StochasticInteraction("NuclearDecay") {
	std::string filename = getDataPath("/NuclearDecay/decay_table.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("mpc::NuclearDecay: could not open file " + filename);

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
	return "Nuclear decay";
}

bool NuclearDecay::setNextInteraction(Candidate *candidate) const {
	int id = candidate->current.getId();

	std::map<int, std::vector<InteractionState> >::const_iterator iMode =
			decayTable.find(id);
	if (iMode == decayTable.end())
		return false;

	const std::vector<InteractionState> &states = iMode->second;

	// find decay mode with minimum random decay distance
	InteractionState decay;
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
	candidate->setInteractionState(name, decay);
	return true;
}

void NuclearDecay::performInteraction(Candidate *candidate) const {
	InteractionState decay;
	candidate->getInteractionState(name, decay);
	candidate->clearInteractionStates();

	// parse decay channel
	int nBeta = digit(decay.channel, 10000);
	int nBetaPlus = digit(decay.channel, 1000);
	int nAlpha = digit(decay.channel, 100);
	int nProton = digit(decay.channel, 10);
	int nNeutron = digit(decay.channel, 1);

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
		candidate->addSecondary(getNucleusId(4, 2), EpA * 4);
	}
	for (size_t i = 0; i < nProton; i++) {
		candidate->addSecondary(getNucleusId(1, 1), EpA);
	}
	for (size_t i = 0; i < nNeutron; i++) {
		candidate->addSecondary(getNucleusId(1, 0), EpA);
	}
}

} // namespace mpc

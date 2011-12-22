#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/Module.h"
#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/MersenneTwister.h"

#include <math.h>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>

namespace mpc {

struct DecayMode {
	int channel;
	double distance;
};

class Decay: public Module {
private:
	MTRand mtrand;
	std::map<int, std::vector<DecayMode> > decayModeMap;
	int cached_id;
	int cached_channel;
	double cached_distance;

public:
	Decay() {
		cached_id = 0;
		// read decay data
		int id;
		DecayMode dm;
		std::ifstream infile("data/Decay/decayTable.txt");
		while (infile.good()) {
			if (infile.peek() != '#') {
				infile >> id >> dm.distance >> dm.channel;
				if (infile) {
					dm.distance *= c_light;
					decayModeMap[id].push_back(dm);
				}
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		infile.close();
	}

	std::string getDescription() const {
		return "Nuclear Decay";
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		int id = candidate->current.getId();
		std::vector<DecayMode> decayModes = decayModeMap[id];
		double gamma = candidate->current.getLorentzFactor();
		double step = candidate->getCurrentStep() / gamma;

		// check if stable
		if (decayModes.size() == 0)
			return;

		// check if a decay is already cached, if not select one
		if (id != cached_id) {
			cached_id = id;
			// find decay mode with minimum random decay distance
			cached_distance = std::numeric_limits<double>::max();
			for (size_t i = 0; i < decayModes.size(); i++) {
				double d = -log(mtrand.rand()) * decayModes[i].distance;
				if (d > cached_distance)
					continue;
				cached_distance = d;
				cached_channel = decayModes[i].channel;
			}
		}

		// check if life-time is over, reduce life-time
		if (cached_distance > step) {
			cached_distance -= step;
			candidate->limitNextStep(cached_distance * gamma);
			return;
		}

		// else life time is over, decay
		cached_id = 0;
		doDecay(candidate, secondaries);
	}

	void doDecay(Candidate *candidate, std::vector<Candidate *> &secondaries) {
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
};

} // namespace mpc

#endif /* DECAY_H_ */

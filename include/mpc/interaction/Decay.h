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
		// read decay data
		std::ifstream infile("include/data/decay.dat");
		char header[256];
		infile.getline(header, 255);

		int id;
		DecayMode dm;
		while (infile.good()) {
			infile >> id >> dm.distance >> dm.channel;
			decayModeMap[id].push_back(dm);
		}
	}

	std::string getDescription() const {
		return "Nuclear Decay";
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		int id = candidate->current.getId();
		double gamma = candidate->current.getLorentzFactor();
		double step = candidate->getCurrentStep() / gamma;
		std::vector<DecayMode> decayModes = decayModeMap[id];

		std::cout << id << ", A " << candidate->current.getMassNumber() << ", Z " << candidate->current.getChargeNumber()<< std::endl;

		// check if stable
		if (decayModes.size() == 0)
			std::cout << "stable" << std::endl;
			return;

		// check if decay is already cached, if not select decay
		if (id != cached_id) {
			id = cached_id;
			// find decay mode with minimum random decay distance
			cached_distance = std::numeric_limits<double>::max();
			for (size_t i = 0; i < decayModes.size(); i++) {
				double d = -log(mtrand.rand()) * decayModes[i].distance;
				std::cout << d / Mpc << std::endl;
				if (d > cached_distance)
					continue;
				cached_distance = d;
				cached_channel = decayModes[i].channel;
			}
		}

		// check if life-time is over
		if (cached_distance > step) {
			cached_distance -= step;
			candidate->limitNextStep(cached_distance * gamma);
			return;
		}

		// life time over -> decay
		int nBeta = cached_channel / 10000;
		int nBetaPlus = (cached_channel % 10000) / 1000;
		int nAlpha = (cached_channel % 1000) / 100;
		int nProton = (cached_channel % 100) / 10;
		int nNeutron = cached_channel % 10;

		// update particle
		int dA = 4 * nAlpha + nProton + nNeutron;
		int dZ = 2 * nAlpha + nProton;

		int A = getMassNumberFromNucleusId(id);
		int Z = getChargeNumberFromNucleusId(id);
		double energyPerNucleon = candidate->current.getEnergy() / double(A);

		candidate->current.setId(getNucleusId(A - dA, Z - dZ));
		candidate->current.setEnergy(energyPerNucleon * (A - dA));

		// create secondaries
		for (int i = 0; i < nBeta; i++) {
			// electron + neutrino
		}
		for (int i = 0; i < nBetaPlus; i++) {
			// positron + neutrino
		}
		for (int i = 0; i < nAlpha; i++) {
			createSecondary(candidate, secondaries, getNucleusId(4, 2),
					energyPerNucleon * 4);
		}
		for (int i = 0; i < nProton; i++) {
			createSecondary(candidate, secondaries, getNucleusId(4, 2),
					energyPerNucleon);
		}
		for (int i = 0; i < nNeutron; i++) {
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

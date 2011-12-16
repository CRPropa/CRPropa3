#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/Module.h"
#include "mpc/MersenneTwister.h"

#include <math.h>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>

namespace mpc {

struct DecayMode {
	int decayType;
	double decayTime;
	double deltaMass;
};

class Decay: public Module {
private:
	MTRand mtrand;
	std::map<size_t, std::vector<DecayMode> > decayModeMap;
	size_t cached_id; // id of particle for which the decay was set
	DecayMode cached_decayMode; // set decay mode
	double cashed_distance; // decay

public:
	Decay() {
		cashed_distance(0), cached_id(0);

		// read decay table
		char header[256];
		std::ifstream infile("include/data/bnl_data.dat");
		infile.getline(header, 255);

		int id;
		DecayMode dm;
		while (!infile.eof()) {
			infile >> id >> dm.decayType >> dm.decayTime;
			decayModeMap[id].push_back(dm);
		}
	}

	std::string getDescription() const {
		return "Nuclear Decay";
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		int id = candidate.current.getId();
		double gamma = candidate.current.getLorentzFactor()
		double step = candidate.getCurrentStep() / gamma;

		std::vector<DecayMode> vecDecayModes = decayModeMap[id];

		if (vecDecayModes.size() == 0) // stable
			return;

		if (id != cached_id) { // new particle, set new decay
			// find decay mode with minimum random decay time
			cashed_distance = std::numeric_limits<double>::max();
			for (size_t i = 0; i < vecDecayModes.size(); i++) {
				double t = -vecDecayModes[i].decayTime * log(mtrand.rand());
				if (t > cashed_distance)
					continue;
				cashed_distance = t;
				cached_decayMode = vecDecayModes[i];
			}
			candidate.limitNextStep(cashed_distance * gamma);
			return;
		}

		if (cashed_distance > step) { // life-time not over
			cashed_distance -= step;
			candidate.limitNextStep(cashed_distance * gamma);
			return;
		}

		// else decay
		switch (cached_decayMode.decayType) {
		case -1:
			std::cout << "Beta- Decay" << std::endl;
			break;
		case 1:
			std::cout << "Beta+ Decay" << std::endl;
			break;
		case 2:
			std::cout << "Proton Dripping" << std::endl;
			break;
		case 3:
			std::cout << "Neutron Dripping" << std::endl;
			break;
		case 4:
			std::cout << "Alpha Decay" << std::endl;
			break;
			// goto middle or recursion?
		}
	}
};

} // namespace mpc

#endif /* DECAY_H_ */

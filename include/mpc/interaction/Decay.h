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
#include <string>
#include <vector>

namespace mpc {

struct DecayMode {
	int decayType;
	double decayTime;
	double deltaMass;
};

class Decay: public Module {
public:
	Decay() {
		initDecayTable();
		timer = 0;
		cached_id = 0;
	}

	void process(Candidate &candidate) {
		size_t id = candidate.current.getId();
		double lf = candidate.current.getLorentzFactor()
		double t = candidate.getCurrentStep() / lf;


		std::vector<DecayMode> vecDecayModes = decayTable[id];
//        getDecayLength(Z,A,tau,daughter,gamma);

		if (vecDecayModes.size() == 0) // stable
			return;

		if (id != cached_id) { // new particle, set new decay
			// find decay mode with minimum random decay time
			timer = std::numeric_limits<double>::max();
			for (size_t i = 0; i < vecDecayModes.size(); i++) {
				double t = -vecDecayModes[i].decayTime * log(mtrand.rand());
				if (t > timer)
					continue;
				timer = t;
				cached_decayMode = vecDecayModes[i];
			}
			// limit next time step
			candidate.setNextStepUpperLimit(timer * lf);
			return;
		}

		if (timer > t) { // life-time not over
			timer -= t;
			// limit next time step
			candidate.setNextStepUpperLimit(timer * lf);
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
			// goto middle
		}
	}

	std::string getDescription() const {
		return "Nuclear Decay";
	}

private:
	MTRand mtrand;
	std::map<size_t, std::vector<DecayMode> > decayTable;
	double timer;
	DecayMode cached_decayMode;
	size_t cached_id;

	void initDecayTable() {
		char str[256];
		size_t A, Z, N;
		std::string name;
		DecayMode dm;

		std::ifstream infile("include/data/bnl_data.dat");
		infile.getline(str, 255);
		infile.getline(str, 255);

		while (!infile.eof()) {
			infile >> A >> Z >> N >> dm.deltaMass >> dm.decayTime
					>> dm.decayType >> name;
			decayTable[size_t(1000000000 + Z * 10000 + A * 10)].push_back(dm);
		}
	}

};

} // namespace mpc

#endif /* DECAY_H_ */

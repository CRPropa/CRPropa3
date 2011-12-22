#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/Module.h"
#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
#include "mpc/MersenneTwister.h"

#include <map>
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
	Decay();
	std::string getDescription() const;
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	void decay(Candidate *candidate, std::vector<Candidate *> &secondaries);
};

} // namespace mpc

#endif /* DECAY_H_ */

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
	int channel; // number of (beta-, beta+, alpha, p, n) decays
	double distance; // decay length in [m]
};

class NuclearDecay: public Module {
private:
	MTRand mtrand;
	std::map<int, std::vector<DecayMode> > modeMap;
	int cached_id;
	int cached_channel;
	double cached_distance;

public:
	NuclearDecay();
	std::string getDescription() const;
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	bool setNextInteraction(Candidate *candidate);
	void performInteraction(Candidate *candidate);
};

} // namespace mpc

#endif /* DECAY_H_ */

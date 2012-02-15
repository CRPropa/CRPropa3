#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/Module.h"
#include "mpc/Random.h"

#include <vector>
#include <map>

namespace mpc {

/**
 @class NuclearDecay
 @brief Nuclear decay of unstable nuclei.

 This module simulates the nuclear decay of unstable nuclei.
 decayTable is a map of nuclei id to (distance, channel) where
 channel is the number of (beta-, beta+, alpha, p, n) decays and
 distance is the mean decay length in [m]
 */
class NuclearDecay: public Module {
private:
	std::string name;
	Random random;
	std::map<int, std::vector<InteractionState> > decayTable;

public:
	NuclearDecay();
	std::string getDescription() const;
	void process(Candidate *candidate);
	bool setNextInteraction(Candidate *candidate);
	void performInteraction(Candidate *candidate);
};

} // namespace mpc

#endif /* DECAY_H_ */

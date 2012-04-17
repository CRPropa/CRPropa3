#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/module/StochasticInteraction.h"
#include "mpc/Random.h"

#include <vector>
#include <map>

namespace mpc {

/**
 @class NuclearDecay
 @brief Nuclear decay of unstable nuclei.

 This module simulates the nuclear decay of unstable nuclei using data from NuDat2.
 */
class NuclearDecay: public StochasticInteraction {
private:
	std::map<int, std::vector<InteractionState> > decayTable;

public:
	NuclearDecay();
	std::string getDescription() const;
	bool setNextInteraction(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;
};

} // namespace mpc

#endif /* DECAY_H_ */

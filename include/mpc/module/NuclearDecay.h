#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/module/StochasticInteraction.h"

#include <vector>

namespace mpc {

/**
 @class NuclearDecay
 @brief Nuclear decay of unstable nuclei.

 This module simulates the nuclear decay of unstable nuclei using data from NuDat2.
 */
class NuclearDecay: public StochasticInteraction {
private:
	bool haveElectrons;
	bool haveNeutrinos;

	std::vector<std::vector<InteractionState> > decayTable;
	// InteractionState.channel is (#beta- #beta+ #alpha #proton #neutron)
	// InteractionState.distance is the mean free path [m]
	double tBeta[50]; // electron kinetic energy [J] in neutron decays
	double cdfBeta[50]; // cumulative distribution function for the electron kinetic energy [J] in neutron decays

public:
	NuclearDecay(bool electrons = false, bool neutrinos = false);
	bool setNextInteraction(Candidate *candidate,
			InteractionState &interaction) const;
	void performInteraction(Candidate *candidate) const;
	void betaDecay(Candidate *candidate, bool isBetaPlus) const;
	void nucleonEmission(Candidate *candidate, int dA, int dZ) const;
};

} // namespace mpc

#endif /* DECAY_H_ */

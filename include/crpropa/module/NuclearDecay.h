#ifndef CRPROPA_DECAY_H
#define CRPROPA_DECAY_H

#include "crpropa/module/StochasticInteraction.h"

#include <vector>

namespace crpropa {

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
	std::vector<double> tBeta; // electron kinetic energy [J] in neutron decays
	std::vector<double> cdfBeta; // cumulative distribution function for the electron kinetic energy [J] in neutron decays

public:
	NuclearDecay(bool electrons = false, bool neutrinos = false);
	void setHaveElectrons(bool b);
	void setHaveNeutrinos(bool b);
	bool setNextInteraction(Candidate *candidate,
			InteractionState &interaction) const;
	void performInteraction(Candidate *candidate) const;
	void betaDecay(Candidate *candidate, bool isBetaPlus) const;
	void nucleonEmission(Candidate *candidate, int dA, int dZ) const;
};

} // namespace crpropa

#endif // CRPROPA_DECAY_H

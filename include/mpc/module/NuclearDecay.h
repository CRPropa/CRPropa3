#ifndef DECAY_H_
#define DECAY_H_

#include "mpc/module/StochasticInteraction.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
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
	gsl_interp_accel *acc;
	gsl_spline *Tbeta; // inverse cdf electron kinetic energy [J] in neutron decay

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

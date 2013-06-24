#ifndef CRPROPA_STOCHASTICINTERACTION_H
#define CRPROPA_STOCHASTICINTERACTION_H

#include "crpropa/Module.h"

namespace crpropa {

/**
 @class StochasticInteraction
 @brief Base class for stochastic interactions
 */
class StochasticInteraction: public Module {
public:
	void process(Candidate *candidate) const;
	virtual bool setNextInteraction(Candidate *candidate,
			InteractionState &interaction) const = 0;
	virtual void performInteraction(Candidate *candidate) const = 0;
};

} // namespace crpropa

#endif // CRPROPA_STOCHASTICINTERACTION_H

#ifndef STOCHASTICINTERACTION_H_
#define STOCHASTICINTERACTION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class StochasticInteraction
 @brief Base class for stochastic interactions
 */
class StochasticInteraction: public Module {
protected:
	std::string name;

public:
	StochasticInteraction(std::string name);
	void process(Candidate *candidate) const;
	virtual bool setNextInteraction(Candidate *candidate) const = 0;
	virtual void performInteraction(Candidate *candidate) const = 0;
};

} // namespace mpc

#endif /* STOCHASTICINTERACTION_H_ */

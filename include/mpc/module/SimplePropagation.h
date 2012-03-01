#ifndef SIMPLEPROPAGATION_H_
#define SIMPLEPROPAGATION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class SimplePropagation
 @brief Simple propagation in absence of magnetic fields.
 */
class SimplePropagation: public Module {
private:
	double acceleration;
	double initialStep;

public:
	SimplePropagation();
	SimplePropagation(double acceleration, double initialStep);
	void process(Candidate *candidate);
	std::string getDescription() const;
};

} // namespace mpc

#endif /* SIMPLEPROPAGATION_H_ */


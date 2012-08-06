#ifndef SIMPLEPROPAGATION_H_
#define SIMPLEPROPAGATION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class SimplePropagation
 @brief Simple rectalinear propagation in absence of magnetic fields.
 */
class SimplePropagation: public Module {
private:
	double acceleration, minStep, maxStep;

public:
	SimplePropagation(double acceleration = 10, double minStep = 0.1 * kpc, double maxStep = 4000 * Mpc);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* SIMPLEPROPAGATION_H_ */


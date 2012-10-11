#ifndef SIMPLEPROPAGATION_H_
#define SIMPLEPROPAGATION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class SimplePropagation
 @brief Simple rectalinear propagation in absence of magnetic fields.

 This module performs a rectalinear propagation.
 It always proposes a next step of the maximum step size.
 */
class SimplePropagation: public Module {
private:
	double minStep, maxStep;

public:
	SimplePropagation(double minStep = 0, double maxStep = 10 * Mpc);
	void process(Candidate *candidate) const;
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);
	std::string getDescription() const;
};

} // namespace mpc

#endif /* SIMPLEPROPAGATION_H_ */


#ifndef SIMPLEPROPAGATION_H_
#define SIMPLEPROPAGATION_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class SimplePropagation
 @brief Simple rectalinear propagation in absence of magnetic fields.

 This module performs a rectalinear propagation.
 The step size is guaranteed to be larger than minStep and smaller than maxStep.
 It always proposes a next step size of maxStep.
 */
class SimplePropagation: public Module {
private:
	double minStep, maxStep;

public:
	SimplePropagation(double minStep = 0, double maxStep = 10 * Mpc);
	void process(Candidate *candidate) const;
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);
	double getMinimumStep() const;
	double getMaximumStep() const;
	std::string getDescription() const;
};

} // namespace mpc

#endif // SIMPLEPROPAGATION_H_


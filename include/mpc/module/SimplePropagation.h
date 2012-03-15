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
	double acceleration;
	double minimumStep;

public:
	SimplePropagation(double acceleration = 10, double minimumStep = 0.1 * kpc);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* SIMPLEPROPAGATION_H_ */


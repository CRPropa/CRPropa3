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
public:
	SimplePropagation();
	SimplePropagation(double acceleration);
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	std::string getDescription() const;
};

} // namespace mpc

#endif /* SIMPLEPROPAGATION_H_ */


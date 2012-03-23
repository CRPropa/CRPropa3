#ifndef MPC_MODULE_H_
#define MPC_MODULE_H_

#include "mpc/Candidate.h"
#include "mpc/Referenced.h"

#include <string>

namespace mpc {

class Candidate;

/**
 @class Module
 @brief Abstract base class for modules
 */
class Module: public Referenced {
public:
	virtual ~Module() {
	}

	virtual std::string getDescription() const;

	virtual void process(Candidate *candidate) const = 0;

	void process(ref_ptr<Candidate> &candidate) const {
		process(candidate.get());
	}
};

} // namespace mpc

#endif /* MPC_MODULE_H_ */

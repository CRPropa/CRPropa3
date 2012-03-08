#ifndef MPC_MODULE_H_
#define MPC_MODULE_H_

#include "mpc/Candidate.h"
#include "mpc/Referenced.h"

#include <string>
#include <vector>

namespace mpc {

class Candidate;

/**
 @class Module
 @brief Module base class.
 */
class Module: public Referenced {
public:
	virtual ~Module() {
	}

	virtual std::string getDescription() const;

	virtual void process(Candidate *candidate) = 0;
};

} // namespace mpc

#endif /* MPC_MODULE_H_ */

#ifndef MODULE_H_
#define MODULE_H_

#include "mpc/Candidate.h"
#include <string>
#include <vector>

namespace mpc {

class Candidate;

/**
 @class Module
 @brief Module base class.
 */
class Module {
public:
	virtual ~Module() {
	}

	virtual std::string getDescription() const;

	virtual void process(Candidate *candidate,
			std::vector<Candidate *> &secondaries) = 0;
};

}

#endif /* MODULE_H_ */

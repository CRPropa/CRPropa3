#ifndef MODULE_H_
#define MODULE_H_

#include <string>
#include <vector>

namespace mpc {

class Candidate;

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

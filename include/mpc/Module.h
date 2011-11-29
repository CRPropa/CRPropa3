#ifndef MODULE_H_
#define MODULE_H_

#include <string>

namespace mpc {

class Candidate;

class Module {
public:
	virtual ~Module() {
	}

	virtual std::string getDescription() const;

	virtual void apply(Candidate &candidate) = 0;
};

}

#endif /* MODULE_H_ */

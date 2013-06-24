#ifndef CRPROPA_MODULE_H
#define CRPROPA_MODULE_H

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"
#include "crpropa/Common.h"

#include <string>

namespace crpropa {

class Candidate;

/**
 @class Module
 @brief Abstract base class for modules
 */
class Module: public Referenced {
	std::string description;
public:
	Module();
	virtual ~Module() {
	}
	virtual std::string getDescription() const;
	void setDescription(const std::string &description);
	virtual void process(Candidate *candidate) const = 0;
	inline void process(ref_ptr<Candidate> candidate) const {
		process(candidate.get());
	}
};

} // namespace crpropa

#endif /* CRPROPA_MODULE_H */

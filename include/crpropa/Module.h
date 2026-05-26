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
class Module {
	std::string description;
public:
	Module();
	virtual ~Module() {
	}
	virtual std::string getDescription() const;
	void setDescription(const std::string &description);
	virtual void process(ref_ptr<Candidate> candidate) const = 0;
};


/**
 @class AbstractCondition
 @brief Abstract Module providing common features for conditional modules.
 */
class AbstractCondition: public Module {
protected:
	ref_ptr<Module> rejectAction, acceptAction;
	bool makeRejectedInactive, makeAcceptedInactive;
	std::string rejectFlagKey, rejectFlagValue;
	std::string acceptFlagKey, acceptFlagValue;

	void reject(ref_ptr<Candidate> candidate) const;
	void accept(ref_ptr<Candidate> candidate) const;

public:
	AbstractCondition();
	void onReject(ref_ptr<Module> rejectAction);
	void onAccept(ref_ptr<Module> acceptAction);
	void setMakeRejectedInactive(bool makeInactive);
	void setMakeAcceptedInactive(bool makeInactive);
	void setRejectFlag(std::string key, std::string value);
	void setAcceptFlag(std::string key, std::string value);

	// return the reject flag (key & value), delimiter is the "&".
	std::string getRejectFlag();

	// return the accept flag (key & value), delimiter is the "&"
	std::string getAcceptFlag();
};

/**
 @class Deactivation
 @brief Direct deactivation of the candidate. Can be used for debuging.
*/
class Deactivation: public AbstractCondition {
	public: 
		void process(ref_ptr<Candidate> cand) const { reject(cand); }
};


} // namespace crpropa

#endif /* CRPROPA_MODULE_H */

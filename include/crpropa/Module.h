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

 This class is the abstract base for any module, that means
 if your module should be made usable over the ModuleList::add function
 it needs to be derived from Module
 */
class Module: public Referenced {
	std::string description; //< Description of any module, can be overriden by setDescription
public:
	/** Constructor
	 * Automatically sets the description to the name of the derived class
	 */
	Module();
	/** Destructor
	 * Virtual destructor that overrides Referenced::~Referenced, does nothing per default.
	 */
	virtual ~Module() override {
	}
	/// Returns description, can be overriden by derived class
	virtual std::string getDescription() const;
	/** Sets description, can be used instead of getDescription in derived class
	 * @param description  Description string
	 */
	void setDescription(const std::string &description);
	/** Process function
	 * This is the main function to work with, it is called during every ModuleList::run step.
	 * Derived functions need to override this overload of process!
	 * @param candidate  Candidate raw pointer
	 */
	virtual void process(Candidate *candidate) const = 0;
	/** Process function
	 * Overload for process function to explicitly make a conversion from ref_ptr to a raw pointer
	 */
	inline void process(ref_ptr<Candidate> candidate) const {
		process(candidate.get());
	}
	#ifndef SWIG
	/// Insertion operator for all Modules
	friend std::ostream& operator<<(std::ostream& os, const Module& obj) {
		os << obj.getDescription();
		return os;
	}
	#endif
};


/**
 @class AbstractCondition
 @brief Abstract Module providing common features for conditional modules.

 This abstract base class is a expansion on the Module class which only proves basic features.
 Over this class it is possible to set modules that should be called on rejection of the particle
 or on acception.
 */
class AbstractCondition: public Module {
protected:
	ref_ptr<Module> rejectAction; /**< Module to call when Candidate is rejected */
	ref_ptr<Module> acceptAction; /**< Module to call when Candidate is accepted */
	bool makeRejectedInactive; /**< If to make Candidate inactive on rejection */
	bool makeAcceptedInactive; /**< If to make Candidate inactive on acception */
	std::string rejectFlagKey; /**< key for the Candidate property that should be set on rejection */
	std::string rejectFlagValue; /**< value for the Candidate property that should be set on rejection */
	std::string acceptFlagKey; /**< key for the Candidate property that should be set on acception */
	std::string acceptFlagValue; /**< value for the Candidate property that should be set on acception */

	/** Function to reject particle
	 * This function rejects the particle, what that means is determined by
	 * the rejectAction and makeRejectedInactive variable.
	 * @param candidate  Candidate raw pointer to reject
	 */
	void reject(Candidate *candidate) const;
	/** Overload for reject(Candidate *candidate)
	 * This function wraps reject(Candidate *candidate) so a reference pointer
	 * can be conversed explicitly
	 */
	inline void reject(ref_ptr<Candidate> candidate) const {
		reject(candidate.get());
	}

	/** Function to accepts particle
	 * This function accepts the particle, what that means is determined by
	 * the acceptAction and makeAcceptedInactive variable.
	 * @param candidate  Candidate to reject
	 */
	void accept(Candidate *candidate) const;
	/** Overload for accept(Candidate *candidate)
	 * This function wraps accept(Candidate *candidate) so a reference pointer
	 * can be conversed explicitly
	 */
	inline void accept(ref_ptr<Candidate> candidate) const {
		accept(candidate.get());
	}

public:
	/** Default constructor
	 * Sets the following parameter:
	 * @var makeRejectedInactive=true
	 * @var makeAcceptedInactive=false
	 * @var rejectFlagKey="Rejected"
	 * @var rejectFlagValue=typeid(*this).name()
	 */
	AbstractCondition();
	/** Sets the module that should be invoked on rejection
	 * @param rejectAction  Module that should be called on rejection
	 */
	void onReject(Module *rejectAction);
	/** Sets the module that should be invoked on acception
	 * @param rejectAction  Module that should be called on acception
	 */
	void onAccept(Module *acceptAction);
	/** Whether to make Candidate inactive on rejection */
	void setMakeRejectedInactive(bool makeInactive);
	/** Whether to make Candidate inactive on acception */
	void setMakeAcceptedInactive(bool makeInactive);
	/** sets Rejection flag
	 * Sets the key and value that should be set when the Candidate is rejected
	 * @param key  property key
	 * @param value  property value
	 */
	void setRejectFlag(std::string key, std::string value);
	/** sets the Acception flag
	 * Sets the key and value that should be set when the Candidate is accepted
	 * @param key  property key
	 * @param value  property value
	 */
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
		void process(Candidate *cand) const { reject(cand); }
};


} // namespace crpropa

#endif /* CRPROPA_MODULE_H */

#include "crpropa/Module.h"

#include <typeinfo>

namespace crpropa {

Module::Module() {
	const std::type_info &info = typeid(*this);
	setDescription(info.name());
}

std::string Module::getDescription() const {
	return description;
}

void Module::setDescription(const std::string &d) {
	description = d;
}

AbstractCondition::AbstractCondition() :
		makeRejectedInactive(true), makeAcceptedInactive(false), rejectFlagKey(
				"Rejected"), rejectFlagValue( typeid(*this).name() ) {
}

void AbstractCondition::reject(ref_ptr<Candidate> candidate) const {
	if (!candidate.valid())
		return;

	if (rejectAction.valid())
		rejectAction->process(candidate);

	if (!rejectFlagKey.empty())
		candidate->setProperty(rejectFlagKey, rejectFlagValue);

	if (makeRejectedInactive)
		candidate->setActive(false);
}

void AbstractCondition::accept(ref_ptr<Candidate> candidate) const {
	if (!candidate.valid())
		return;

	if (acceptAction.valid())
		acceptAction->process(candidate);

	if (!acceptFlagKey.empty())
		candidate->setProperty(acceptFlagKey, acceptFlagValue);

	if (makeAcceptedInactive)
		candidate->setActive(false);
}

void AbstractCondition::setMakeRejectedInactive(bool deactivate) {
	makeRejectedInactive = deactivate;
}

void AbstractCondition::setMakeAcceptedInactive(bool deactivate) {
	makeAcceptedInactive = deactivate;
}

void AbstractCondition::onReject(ref_ptr<Module> action) {
	rejectAction = action;
}

void AbstractCondition::onAccept(ref_ptr<Module> action) {
	acceptAction = action;
}

void AbstractCondition::setRejectFlag(std::string key, std::string value) {
	rejectFlagKey = key;
	rejectFlagValue = value;
}

void AbstractCondition::setAcceptFlag(std::string key, std::string value) {
	acceptFlagKey = key;
	acceptFlagValue = value;
}

std::string AbstractCondition::getRejectFlag() {
	std::string out = rejectFlagKey + "&" + rejectFlagValue; 
	return out;
}

std::string AbstractCondition::getAcceptFlag() {
	std::string out = acceptFlagKey + "&" + acceptFlagValue;
	return out;
}


} // namespace crpropa

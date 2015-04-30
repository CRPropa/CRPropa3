#include "crpropa/module/BreakCondition.h"
#include "crpropa/Units.h"

#include <sstream>

namespace crpropa {

AbstractBreakCondition::AbstractBreakCondition() :
		makeInactive(true), flagKey("Deactivated") {

}

void AbstractBreakCondition::processBreak(Candidate *candidate) const {
	if (!candidate)
		return;

	if (breakAction.valid())
		breakAction->process(candidate);

	if (!flagKey.empty())
		candidate->setProperty(flagKey, flagValue);

	if (makeInactive)
		candidate->setActive(false);
}

void AbstractBreakCondition::setMakeInactive(bool deactivate) {
	makeInactive = deactivate;
}

void AbstractBreakCondition::onBreak(Module *action) {
	breakAction = action;
}

void AbstractBreakCondition::setFlag(std::string key, std::string value) {
	flagKey = key;
	flagValue = value;
}

void AbstractBreakCondition::endRun() {
	if (breakAction.valid())
		breakAction->endRun();
}

MaximumTrajectoryLength::MaximumTrajectoryLength(double maxLength) :
		maxLength(maxLength) {
}

void MaximumTrajectoryLength::setMaximumTrajectoryLength(double length) {
	maxLength = length;
}

double MaximumTrajectoryLength::getMaximumTrajectoryLength() const {
	return maxLength;
}

std::string MaximumTrajectoryLength::getDescription() const {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc, ";
	s << "Flag: '" << flagKey << "' -> '" << flagValue << "', ";
	s << "MakeInactive: " << (makeInactive ? "yes" : "no");
	if (breakAction.valid())
		s << ", Action: " << breakAction->getDescription();
	return s.str();
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double l = c->getTrajectoryLength();
	if (l >= maxLength) {
		processBreak(c);
	} else {
		c->limitNextStep(maxLength - l);
	}
}

MinimumEnergy::MinimumEnergy(double minEnergy) :
		minEnergy(minEnergy) {
}

void MinimumEnergy::setMinimumEnergy(double energy) {
	minEnergy = energy;
}

double MinimumEnergy::getMinimumEnergy() const {
	return minEnergy;
}

void MinimumEnergy::process(Candidate *c) const {
	if (c->current.getEnergy() > minEnergy)
		return;
	else
		processBreak(c);
}

std::string MinimumEnergy::getDescription() const {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV, ";
	s << "Flag: '" << flagKey << "' -> '" << flagValue << "', ";
	s << "MakeInactive: " << (makeInactive ? "yes" : "no");
	if (breakAction.valid())
		s << ", Action: " << breakAction->getDescription();
	return s.str();
}

MinimumRedshift::MinimumRedshift(double zmin) :
		zmin(zmin) {
}

void MinimumRedshift::setMinimumRedshift(double z) {
	zmin = z;
}

double MinimumRedshift::getMinimumRedshift() {
	return zmin;
}

void MinimumRedshift::process(Candidate* c) const {
	if (c->getRedshift() > zmin)
		return;
	else
		processBreak(c);
}

std::string MinimumRedshift::getDescription() const {
	std::stringstream s;
	s << "Minimum redshift: " << zmin << ", ";
	s << "Flag: '" << flagKey << "' -> '" << flagValue << "', ";
	s << "Meactivate: " << makeInactive ? "yes" : "no";
	if (breakAction.valid())
		s << ", Action: " << breakAction->getDescription();
	return s.str();
}

} // namespace crpropa

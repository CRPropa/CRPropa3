#include "crpropa/module/BreakCondition.h"
#include "crpropa/Units.h"

#include <sstream>

namespace crpropa {

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
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double l = c->getTrajectoryLength();
	if (l >= maxLength) {
		reject(c);
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
		reject(c);
}

std::string MinimumEnergy::getDescription() const {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
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
		reject(c);
}

std::string MinimumRedshift::getDescription() const {
	std::stringstream s;
	s << "Minimum redshift: " << zmin << ", ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

} // namespace crpropa

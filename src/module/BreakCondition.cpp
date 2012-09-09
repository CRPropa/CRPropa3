#include "mpc/module/BreakCondition.h"

#include <sstream>

namespace mpc {

MaximumTrajectoryLength::MaximumTrajectoryLength(double maxLength) {
	this->maxLength = maxLength;
}

void MaximumTrajectoryLength::setMaximumTrajectoryLength(double length) {
	maxLength = length;
}

double MaximumTrajectoryLength::getMaximumTrajectoryLength() const {
	return maxLength;
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double l = c->getTrajectoryLength();
	if (l >= maxLength) {
		c->setActive(false);
		c->setProperty("Deactivated", getDescription());
	} else {
		c->limitNextStep(maxLength - l);
	}
}

std::string MaximumTrajectoryLength::getDescription() const {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc";
	return s.str();
}

MinimumEnergy::MinimumEnergy(double minEnergy) {
	this->minEnergy = minEnergy;
}

void MinimumEnergy::setMinimumEnergy(double energy) {
	minEnergy = energy;
}

double MinimumEnergy::getMinimumEnergy() const {
	return minEnergy;
}

void MinimumEnergy::process(Candidate *c) const {
	if (c->current.getEnergy() <= minEnergy) {
		c->setActive(false);
		c->setProperty("Deactivated", getDescription());
	}
}

std::string MinimumEnergy::getDescription() const {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV";
	return s.str();
}

MinimumRedshift::MinimumRedshift(double z) :
		zmin(z) {
}

void MinimumRedshift::setMinimumRedshift(double z) {
	zmin = z;
}

double MinimumRedshift::getMinimumRedshift() {
	return zmin;
}

void MinimumRedshift::process(Candidate* c) const {
	double z = c->getRedshift();
	if (z <= zmin) {
		c->setActive(false);
		c->setProperty("Deactivated", getDescription());
	}
}

std::string MinimumRedshift::getDescription() const {
	std::stringstream s;
	s << "Minimum redshift: " << zmin;
	return s.str();
}

} // namespace mpc

#include "mpc/module/BreakCondition.h"

#include <sstream>

namespace mpc {

MaximumTrajectoryLength::MaximumTrajectoryLength(double maxLength) {
	this->maxLength = maxLength;
	updateDescription();
}

void MaximumTrajectoryLength::setMaximumTrajectoryLength(double length) {
	maxLength = length;
	updateDescription();
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

void MaximumTrajectoryLength::updateDescription() {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc";
	setDescription(s.str());
}

MinimumEnergy::MinimumEnergy(double minEnergy) {
	this->minEnergy = minEnergy;
	updateDescription();
}

void MinimumEnergy::setMinimumEnergy(double energy) {
	minEnergy = energy;
	updateDescription();
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

void MinimumEnergy::updateDescription() {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV";
	setDescription(s.str());
}

} // namespace mpc

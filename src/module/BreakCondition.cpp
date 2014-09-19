#include "crpropa/module/BreakCondition.h"
#include "crpropa/Units.h"

#include <sstream>

namespace crpropa {

MaximumTrajectoryLength::MaximumTrajectoryLength(double maxLength,
		std::string flag) :
		maxLength(maxLength), flag(flag) {
}

void MaximumTrajectoryLength::setMaximumTrajectoryLength(double length) {
	maxLength = length;
}

double MaximumTrajectoryLength::getMaximumTrajectoryLength() const {
	return maxLength;
}

void MaximumTrajectoryLength::setFlag(std::string f) {
	flag = f;
}

std::string MaximumTrajectoryLength::getFlag() const {
	return flag;
}

std::string MaximumTrajectoryLength::getDescription() const {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc, flag: "
			<< flag;
	return s.str();
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double l = c->getTrajectoryLength();
	if (l >= maxLength) {
		c->setActive(false);
		c->setProperty(flag, getDescription());
	} else {
		c->limitNextStep(maxLength - l);
	}
}

MinimumEnergy::MinimumEnergy(double minEnergy, std::string flag) :
		minEnergy(minEnergy), flag(flag) {
}

void MinimumEnergy::setMinimumEnergy(double energy) {
	minEnergy = energy;
}

double MinimumEnergy::getMinimumEnergy() const {
	return minEnergy;
}

void MinimumEnergy::setFlag(std::string f) {
	flag = f;
}

std::string MinimumEnergy::getFlag() const {
	return flag;
}

void MinimumEnergy::process(Candidate *c) const {
	if (c->current.getEnergy() > minEnergy)
		return;
	c->setActive(false);
	c->setProperty(flag, getDescription());
}

std::string MinimumEnergy::getDescription() const {
	std::stringstream s;
	s << "Minimum energy: " << minEnergy / EeV << " EeV, flag: " << flag;
	return s.str();
}

MinimumRedshift::MinimumRedshift(double zmin, std::string flag) :
		zmin(zmin), flag(flag) {
}

void MinimumRedshift::setMinimumRedshift(double z) {
	zmin = z;
}

double MinimumRedshift::getMinimumRedshift() {
	return zmin;
}

void MinimumRedshift::setFlag(std::string f) {
	flag = f;
}

std::string MinimumRedshift::getFlag() const {
	return flag;
}

void MinimumRedshift::process(Candidate* c) const {
	if (c->getRedshift() > zmin)
		return;
	c->setActive(false);
	c->setProperty(flag, getDescription());
}

std::string MinimumRedshift::getDescription() const {
	std::stringstream s;
	s << "Minimum redshift: " << zmin << ", flag: " << flag;
	return s.str();
}

} // namespace crpropa

#include "crpropa/module/BreakCondition.h"
#include "crpropa/ParticleID.h"
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

void MaximumTrajectoryLength::addObserverPosition(const Vector3d& position) {
	observerPositions.push_back(position);
}

const std::vector<Vector3d>& MaximumTrajectoryLength::getObserverPositions() const {
	return observerPositions;
}

std::string MaximumTrajectoryLength::getDescription() const {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / Mpc << " Mpc, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	s << "\n  Observer positions: \n";
	for (size_t i = 0; i < observerPositions.size(); i++)
		s << "    - " << observerPositions[i] / Mpc << " Mpc\n";
	return s.str();
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double length = c->getTrajectoryLength();
	Vector3d position = c->current.getPosition();

	if(observerPositions.size()) {
		bool inRange = false;
		for (size_t i = 0; i < observerPositions.size(); i++) {
			double distance = position.getDistanceTo(observerPositions[i]);
			if (distance + length < maxLength)
				inRange = true;
		}
		if (!inRange) {
			reject(c);
			return;
		}
	}

	if (length >= maxLength) {
		reject(c);
	} else {
		c->limitNextStep(maxLength - length);
	}
}

//*****************************************************************************
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

//*****************************************************************************
MinimumRigidity::MinimumRigidity(double minRigidity) :
		minRigidity(minRigidity) {
}

void MinimumRigidity::setMinimumRigidity(double minRigidity) {
	this->minRigidity = minRigidity;
}

double MinimumRigidity::getMinimumRigidity() const {
	return minRigidity;
}

void MinimumRigidity::process(Candidate *c) const {
	double rigidity = fabs(c->current.getEnergy() / chargeNumber(c->current.getId()));
	if (rigidity < minRigidity)
		reject(c);
}

std::string MinimumRigidity::getDescription() const {
	std::stringstream s;
	s << "Minimum rigidity: " << minRigidity / EeV << " EeV, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

//*****************************************************************************
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

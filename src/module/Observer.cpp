#include "crpropa/module/Observer.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Cosmology.h"

#include "kiss/logger.h"

#include <iostream>
#include <cmath>

namespace crpropa {

// Observer -------------------------------------------------------------------
Observer::Observer() :
		makeInactive(true), clone(false) {
}

void Observer::add(ObserverFeature *feature) {
	features.push_back(feature);
}

void Observer::onDetection(Module *action, bool clone_) {
	detectionAction = action;
	clone = clone_;
}

void Observer::process(Candidate *candidate) const {
	// loop over all features and have them check the particle
	DetectionState state = NOTHING;
	for (int i = 0; i < features.size(); i++) {
		DetectionState s = features[i]->checkDetection(candidate);
		if (s == VETO)
			state = VETO;
		else if ((s == DETECTED) && (state != VETO))
			state = DETECTED;
	}

	if (state == DETECTED) {
		for (int i = 0; i < features.size(); i++) {
			features[i]->onDetection(candidate);
		}

		if (detectionAction.valid()) {
			if (clone)
				detectionAction->process(candidate->clone(false));
			else
				detectionAction->process(candidate);
		}

		if (!flagKey.empty())
			candidate->setProperty(flagKey, flagValue);

		if (makeInactive)
			candidate->setActive(false);
	}
}

void Observer::setFlag(std::string key, std::string value) {
	flagKey = key;
	flagValue = value;
}

std::string Observer::getDescription() const {
	std::stringstream ss;
	ss << "Observer";
	for (int i = 0; i < features.size(); i++)
		ss << "\n    " << features[i]->getDescription() << "\n";
	ss << "    Flag: '" << flagKey << "' -> '" << flagValue << "'\n";
	ss << "    MakeInactive: " << (makeInactive ? "yes\n" : "no\n");
	if (detectionAction.valid())
		ss << "    Action: " << detectionAction->getDescription() << ", clone: " << (clone ? "yes" : "no");

	return ss.str();
}

void Observer::setDeactivateOnDetection(bool deactivate) {
	makeInactive = deactivate;
}

// ObserverFeature ------------------------------------------------------------
DetectionState ObserverFeature::checkDetection(Candidate *candidate) const {
	return NOTHING;
}

void ObserverFeature::onDetection(Candidate *candidate) const {
}

std::string ObserverFeature::getDescription() const {
	return description;
}

// ObserverDetectAll ----------------------------------------------------------
DetectionState ObserverDetectAll::checkDetection(Candidate *candidate) const {
	return DETECTED;
}

std::string ObserverDetectAll::getDescription() const {
	return description;
}

// ObserverTracking --------------------------------------------------------
ObserverTracking::ObserverTracking(Vector3d center, double radius, double stepSize) :
		center(center), radius(radius), stepSize(stepSize) {
	if (stepSize == 0) {
		stepSize = radius / 10.;
	}
}

DetectionState ObserverTracking::checkDetection(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getR();

	// no detection if outside of observer sphere
	if (d > radius) {
		// conservatively limit next step to prevent overshooting
		candidate->limitNextStep(fabs(d - radius));

		return NOTHING;
	} else {
		// limit next step
		candidate->limitNextStep(stepSize);

		return DETECTED;
	}
}

std::string ObserverTracking::getDescription() const {
	std::stringstream ss;
	ss << "ObserverTracking: ";
	ss << "center = " << center / Mpc << " Mpc, ";
	ss << "radius = " << radius / Mpc << " Mpc";
	ss << "stepSize = " << stepSize / Mpc << " Mpc";
	return ss.str();
}

// ObserverPoint --------------------------------------------------------------
DetectionState ObserverPoint::checkDetection(Candidate *candidate) const {
	KISS_LOG_WARNING << "ObserverPoint is deprecated and is no longer supported. Please use Observer1D instead.\n";
	double x = candidate->current.getPosition().x;
	if (x > 0) {
		// Limits the next step size to prevent candidates from overshooting in case of non-detection
		candidate->limitNextStep(x);
		return NOTHING;
	}
	// Detects particles when reaching x = 0
	return DETECTED;
}

std::string ObserverPoint::getDescription() const {
	return "ObserverPoint: observer at x = 0";
}

// Observer1D --------------------------------------------------------------
DetectionState Observer1D::checkDetection(Candidate *candidate) const {
	double x = candidate->current.getPosition().x;
	if (x > 0) {
		// Limits the next step size to prevent candidates from overshooting in case of non-detection
		candidate->limitNextStep(x);
		return NOTHING;
	}
	// Detects particles when reaching x = 0
	return DETECTED;
}

std::string Observer1D::getDescription() const {
	return "Observer1D: observer at x = 0";
}

// ObserverRedshiftWindow -----------------------------------------------------
ObserverRedshiftWindow::ObserverRedshiftWindow(double zmin, double zmax) :
		zmin(zmin), zmax(zmax) {
}

DetectionState ObserverRedshiftWindow::checkDetection(
		Candidate *candidate) const {
	double z = candidate->getRedshift();
	if (z > zmax)
		return VETO;
	if (z < zmin)
		return VETO;
	return NOTHING;
}

std::string ObserverRedshiftWindow::getDescription() const {
	std::stringstream ss;
	ss << "ObserverRedshiftWindow: z = " << zmin << " - " << zmax;
	return ss.str();
}

// ObserverInactiveVeto -------------------------------------------------------
DetectionState ObserverInactiveVeto::checkDetection(Candidate *c) const {
	if (not(c->isActive()))
		return VETO;
	return NOTHING;
}

std::string ObserverInactiveVeto::getDescription() const {
	return "ObserverInactiveVeto";
}

// ObserverNucleusVeto --------------------------------------------------------
DetectionState ObserverNucleusVeto::checkDetection(Candidate *c) const {
	if (isNucleus(c->current.getId()))
		return VETO;
	return NOTHING;
}

std::string ObserverNucleusVeto::getDescription() const {
	return "ObserverNucleusVeto";
}

// ObserverNeutrinoVeto -------------------------------------------------------
DetectionState ObserverNeutrinoVeto::checkDetection(Candidate *c) const {
	int id = abs(c->current.getId());
	if ((id == 12) or (id == 14) or (id == 16))
		return VETO;
	return NOTHING;
}

std::string ObserverNeutrinoVeto::getDescription() const {
	return "ObserverNeutrinoVeto";
}

// ObserverPhotonVeto ---------------------------------------------------------
DetectionState ObserverPhotonVeto::checkDetection(Candidate *c) const {
	if (c->current.getId() == 22)
		return VETO;
	return NOTHING;
}

std::string ObserverPhotonVeto::getDescription() const {
	return "ObserverPhotonVeto";
}

// ObserverElectronVeto ---------------------------------------------------------
DetectionState ObserverElectronVeto::checkDetection(Candidate *c) const {
	if (abs(c->current.getId()) == 11)
		return VETO;
	return NOTHING;
}

std::string ObserverElectronVeto::getDescription() const {
	return "ObserverElectronVeto";
}

// ObserverCustomVeto -------------------------------------------------------
ObserverParticleIdVeto::ObserverParticleIdVeto(int pId) {
	vetoParticleId = pId;
}

DetectionState ObserverParticleIdVeto::checkDetection(Candidate *c) const {
	if (c->current.getId() == vetoParticleId)
		return VETO;
	return NOTHING;
}

std::string ObserverParticleIdVeto::getDescription() const {
	return "ObserverParticleIdVeto";
}


// ObserverTimeEvolution --------------------------------------------------------
ObserverTimeEvolution::ObserverTimeEvolution() {}

ObserverTimeEvolution::ObserverTimeEvolution(double min, double dist, double numb) {
	double max = min + numb * dist;
	bool log = false;
	addTimeRange(min, max, numb, log);
}

ObserverTimeEvolution::ObserverTimeEvolution(double min, double max, double numb, bool log) {
	addTimeRange(min, max, numb, log);
}


DetectionState ObserverTimeEvolution::checkDetection(Candidate *c) const {

	if (detList.size()) {
		double length = c->getTrajectoryLength();
		size_t index;
		const std::string DI = "DetectionIndex";
		std::string value;

		// Load the last detection index
		if (c->hasProperty(DI)) {
			index = c->getProperty(DI).asUInt64();
		}
		else {
			index = 0;
		}

		// Break if the particle has been detected once for all detList entries.
		if (index > detList.size()) {
			return NOTHING;
		}

		// Calculate the distance to next detection
		double distance = length - detList[index];

		// Limit next step and detect candidate.
		// Increase the index by one in case of detection
		if (distance < 0.) {
			c->limitNextStep(-distance);
			return NOTHING;
		}
		else {

			if (index < detList.size()-1) {
				c->limitNextStep(detList[index+1]-length);
			}
			c->setProperty(DI, Variant::fromUInt64(index+1));

			return DETECTED;
		}
	}
	return NOTHING;
}

void ObserverTimeEvolution::addTime(const double& t) {
	detList.push_back(t);
}

void ObserverTimeEvolution::addTimeRange(double min, double max, double numb, bool log) {
	for (size_t i = 0; i < numb; i++) {
		if (log == true) {
			addTime(min * pow(max / min, i / (numb - 1.0)));
		} else {
			addTime(min + i * (max - min) / numb);
		}
	}
}

const std::vector<double>& ObserverTimeEvolution::getTimes() const {
	return detList;
}

std::string ObserverTimeEvolution::getDescription() const {
	std::stringstream s;
	s << "List of Detection lengths in kpc";
	for (size_t i = 0; i < detList.size(); i++)
	  s << "  - " << detList[i] / kpc;
	return s.str();
}

// ObserverSurface--------------------------------------------------------------
ObserverSurface::ObserverSurface(Surface* _surface) : surface(_surface) { }

DetectionState ObserverSurface::checkDetection(Candidate *candidate) const
{
		double currentDistance = surface->distance(candidate->current.getPosition());
		double previousDistance = surface->distance(candidate->previous.getPosition());
		candidate->limitNextStep(fabs(currentDistance));

		if (currentDistance * previousDistance > 0)
			return NOTHING;
		else if (previousDistance == 0)
			return NOTHING;
		else
			return DETECTED;
}

std::string ObserverSurface::getDescription() const {
	std::stringstream ss;
	ss << "ObserverSurface: << " << surface->getDescription();
	return ss.str();
}

} // namespace crpropa

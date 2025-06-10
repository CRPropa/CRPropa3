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
	setIsLogarithmicScaling(false);
	setMaximum(min + (numb - 1) * dist);
	setMinimum(min);
	setNIntervals(numb);
}

ObserverTimeEvolution::ObserverTimeEvolution(double min, double max, double numb, bool log) {
	setIsLogarithmicScaling(log);
	setMinimum(min);
	setMaximum(max);
	setNIntervals(numb);
}

ObserverTimeEvolution::ObserverTimeEvolution(const std::vector<double> &detList){
	setTimes(detList);
}

DetectionState ObserverTimeEvolution::checkDetection(Candidate *c) const {

	if (nIntervals) {
		double length = c->getTrajectoryLength();
		size_t index;
		const std::string DI = "DetectionIndex";

		// Load the last detection index
		if (c->hasProperty(DI)) {
			index = c->getProperty(DI).asUInt64();
		}
		else {
			index = 0;
		}

		// Break if the particle has been detected once for all possible times.
		if (index >= nIntervals) {
			return NOTHING;
		}

		// Calculate the distance to next detection
		double distance = length - getTime(index);

		// Limit next step and detect candidate.
		// Increase the index by one in case of detection
		if (distance < 0.) {
			c->limitNextStep(-distance);
			return NOTHING;
		}
		else {

			if (index < nIntervals-2) {
				c->limitNextStep(getTime(index+1)-length);
			}
			c->setProperty(DI, Variant::fromUInt64(index+1));

			return DETECTED;
		}
	}
	return NOTHING;
}

void ObserverTimeEvolution::clear(){
	doDetListConstruction = false;
	detList.clear();
	detList.resize(0);
	setNIntervals(0);
}

void ObserverTimeEvolution::constructDetListIfEmpty(){
	if (detList.empty() && doDetListConstruction) {
		std::vector<double> detListTemp;
		size_t counter = 0;
		while (getTime(counter)<=maximum) {
			detListTemp.push_back(getTime(counter));
			counter++;
		}
		detList.assign(detListTemp.begin(), detListTemp.end());
	}
}

void ObserverTimeEvolution::addTime(const double &time){
	constructDetListIfEmpty();
	detList.push_back(time);
	setNIntervals(nIntervals + 1);  // increase number of entries by one
}

void ObserverTimeEvolution::addTimeRange(double min, double max, double numb, bool log) {
	for (size_t i = 0; i < numb; i++) {
		if (log) {
			if ( min <= 0 ){
				std::cout << "min can not be <= 0 if log=true" << std::endl;
				throw new std::runtime_error("min can not be <= 0 if log=true");
			}
			addTime(min * pow(max / min, i / (numb - 1.0)));
		} else {
			addTime(min + i * (max - min) / (numb - 1.0));
		}
	}
	// allready corrected by addTime, just here to be safe
	setNIntervals(detList.size());
}

void ObserverTimeEvolution::setTimes(const std::vector<double> &detList){
	this->detList.assign(detList.begin(), detList.end());
	setNIntervals(detList.size());
	setMinimum(detList.front());
	setMaximum(detList.back());
	doDetListConstruction = false;
}

void ObserverTimeEvolution::setMinimum(double min){
	if ( (min <= 0) && isLogarithmicScaling){
		std::cout << "minimum can not be <= 0 if isLogarithmicScaling=true" << std::endl;
		throw new std::runtime_error("minimum can not be <= 0 if isLogarithmicScaling=true");
	}
	this->minimum = min;
}

double ObserverTimeEvolution::getTime(size_t index) const {
	if (!detList.empty()) {
		return detList.at(index);
	} else if (isLogarithmicScaling) {
		return minimum * pow(maximum / minimum, index / (nIntervals - 1.0));
	} else {
		return minimum + index * (maximum - minimum) / (nIntervals - 1.0);
	}
}

const std::vector<double>& ObserverTimeEvolution::getTimes() const {
	tempDetList.resize(nIntervals);
	for (size_t i = 0; i < nIntervals; i++){
		tempDetList[i] = getTime(i);
	}
	return tempDetList;
}

std::string ObserverTimeEvolution::getDescription() const {
	std::stringstream s;
	s << "List of Detection lengths in kpc";
	for (size_t i = 0; i < nIntervals; i++)
	  s << "  - " << getTime(i) / kpc;
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

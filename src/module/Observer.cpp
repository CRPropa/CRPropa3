#include "crpropa/module/Observer.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Cosmology.h"

namespace crpropa {

// Observer -------------------------------------------------------------------
Observer::Observer() :
		makeInactive(true), clone(false) {
}

void Observer::add(ObserverFeature *feature) {
	features.push_back(feature);
}

void Observer::beginRun() {
	for (int i = 0; i < features.size(); i++)
		features[i]->beginRun();
}

void Observer::endRun() {
	for (int i = 0; i < features.size(); i++)
		features[i]->endRun();

	if (detectionAction.valid())
		detectionAction->endRun();
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

void ObserverFeature::beginRun() {
}

void ObserverFeature::endRun() {
}

std::string ObserverFeature::getDescription() const {
	return description;
}

// ObserverSmallSphere --------------------------------------------------------
ObserverSmallSphere::ObserverSmallSphere(Vector3d center, double radius) :
		center(center), radius(radius) {
}

DetectionState ObserverSmallSphere::checkDetection(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getR();

	// conservatively limit next step to prevent overshooting
	candidate->limitNextStep(fabs(d - radius));

	// no detection if outside of observer sphere
	if (d > radius)
		return NOTHING;

	// previous distance to observer sphere center
	double dprev = (candidate->previous.getPosition() - center).getR();

	// if particle was inside of sphere in previous step it has already been detected
	if (dprev <= radius)
		return NOTHING;

	// else detection
	return DETECTED;
}

std::string ObserverSmallSphere::getDescription() const {
	std::stringstream ss;
	ss << "ObserverSmallSphere: ";
	ss << "center = " << center / Mpc << " Mpc, ";
	ss << "radius = " << radius / Mpc << " Mpc";
	return ss.str();
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

// ObserverLargeSphere --------------------------------------------------------
ObserverLargeSphere::ObserverLargeSphere(Vector3d center, double radius) :
		center(center), radius(radius) {
}

DetectionState ObserverLargeSphere::checkDetection(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getR();

	// conservatively limit next step size to prevent overshooting
	candidate->limitNextStep(fabs(radius - d));

	// no detection if inside observer sphere
	if (d < radius)
		return NOTHING;

	// previous distance to observer sphere center
	double dprev = (candidate->previous.getPosition() - center).getR();

	// if particle was outside of sphere in previous step it has already been detected
	if (dprev >= radius)
		return NOTHING;

	// else: detection
	return DETECTED;
}

std::string ObserverLargeSphere::getDescription() const {
	std::stringstream ss;
	ss << "ObserverLargeSphere: ";
	ss << "center = " << center / Mpc << " Mpc, ";
	ss << "radius = " << radius / Mpc << " Mpc";
	return ss.str();
}

// ObserverPoint --------------------------------------------------------------
DetectionState ObserverPoint::checkDetection(Candidate *candidate) const {
	double x = candidate->current.getPosition().x;
	if (x > 0) {
		candidate->limitNextStep(x);
		return NOTHING;
	}
	return DETECTED;
}

std::string ObserverPoint::getDescription() const {
	return "ObserverPoint: observer at x = 0";
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



////////////////////////////////////////////////////////////////////////////////

SmallObserverSphere::SmallObserverSphere(Vector3d center, double radius,
		std::string flag, std::string flagValue, bool makeInactive) :
		center(center), radius(radius), flag(flag), flagValue(flagValue), makeInactive(
				makeInactive), maximumTrajectory(
				std::numeric_limits<double>::max()) {
}

void SmallObserverSphere::process(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getR();
	double remainingToBorder = d - radius;

	// conservatively limit next step to prevent overshooting
	candidate->limitNextStep(fabs(remainingToBorder));

	// no detection if outside of observer sphere
	if (d > radius) {

		// disable candidate if it cannot reach this observer
		double remainingToMaximum = maximumTrajectory
				- candidate->getTrajectoryLength();
		if (remainingToBorder > remainingToMaximum)
			candidate->setActive(false);

		return;
	}

	// previous distance to observer sphere center
	double dprev = (candidate->previous.getPosition() - center).getR();

	// if particle was inside of sphere in previous step it has already been detected
	if (dprev <= radius)
		return;

	// else: detection
	candidate->setProperty(flag, flagValue);
	if (makeInactive)
		candidate->setActive(false);
}

void SmallObserverSphere::setCenter(Vector3d c) {
	center = c;
}

void SmallObserverSphere::setRadius(double r) {
	radius = r;
}

void SmallObserverSphere::setFlag(std::string f, std::string v) {
	flag = f;
	flagValue = v;
}

void SmallObserverSphere::setMakeInactive(bool b) {
	makeInactive = b;
}

void SmallObserverSphere::setMaximumTrajectory(double maximumTrajectory) {
	this->maximumTrajectory = maximumTrajectory;
}

std::string SmallObserverSphere::getDescription() const {
	std::stringstream s;
	s << "Small observer sphere: " << radius / Mpc;
	s << " Mpc radius around " << center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	if (maximumTrajectory != std::numeric_limits<double>::max())
		s << ", Maximum trajectory: " << maximumTrajectory / Mpc << "Mpc";
	if (makeInactive)
		s << ", render inactivate";
	return s.str();
}

LargeObserverSphere::LargeObserverSphere(Vector3d center, double radius,
		std::string flag, std::string flagValue, bool makeInactive) :
		center(center), radius(radius), flag(flag), flagValue(flagValue), makeInactive(
				makeInactive) {
}

void LargeObserverSphere::process(Candidate *candidate) const {
	// current distance to observer sphere center
	double d = (candidate->current.getPosition() - center).getR();

	// conservatively limit next step size to prevent overshooting
	candidate->limitNextStep(fabs(radius - d));

	// no detection if inside observer sphere
	if (d < radius)
		return;

	// previous distance to observer sphere center
	double dprev = (candidate->previous.getPosition() - center).getR();

	// if particle was outside of sphere in previous step it has already been detected
	if (dprev >= radius)
		return;

	// else: detection
	candidate->setProperty(flag, flagValue);
	if (makeInactive)
		candidate->setActive(false);
}

void LargeObserverSphere::setCenter(Vector3d c) {
	center = c;
}

void LargeObserverSphere::setRadius(double r) {
	radius = r;
}

void LargeObserverSphere::setFlag(std::string f, std::string v) {
	flag = f;
	flagValue = v;
}

void LargeObserverSphere::setMakeInactive(bool b) {
	makeInactive = b;
}

std::string LargeObserverSphere::getDescription() const {
	std::stringstream s;
	s << "Large observer sphere: " << radius / Mpc;
	s << " Mpc radius around " << center / Mpc;
	s << " Mpc, Flag: '" << flag << "' -> '" << flagValue << "'";
	if (makeInactive)
		s << ", render inactivates";
	return s.str();
}

Observer1D::Observer1D() {
	setDescription("1D observer");
}

void Observer1D::process(Candidate *candidate) const {
	double x = candidate->current.getPosition().x;
	// check if position x > 0
	if (x > std::numeric_limits<double>::min()) {
		candidate->limitNextStep(x);
		return;
	}
	// else: detection
	candidate->setProperty("Detected", "");
	candidate->setActive(false);
}

DetectAll::DetectAll(std::string f, std::string v, bool m) :
		flag(f), flagValue(v), makeInactive(m) {
}

void DetectAll::process(Candidate *candidate) const {
	candidate->setProperty(flag, flagValue);
	if (makeInactive)
		candidate->setActive(false);
}

std::string DetectAll::getDescription() const {
	std::stringstream s;
	s << "DetectAll: Flag: " << flag << " -> " << flagValue;
	if (makeInactive)
		s << ", render inactive";
	return s.str();
}

} // namespace crpropa

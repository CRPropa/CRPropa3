#include "crpropa/module/Observer.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Cosmology.h"

namespace crpropa {

DetectionState ObserverFeature::checkDetection(Candidate *candidate) const {
	return NOTHING;
}

void ObserverFeature::onDetection(Candidate *candidate) const {
}

Observer::Observer(bool makeInactive) :
		makeInactive(makeInactive) {
}

void Observer::add(ObserverFeature *feature) {
	features.push_back(feature);
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

		if (makeInactive)
			candidate->setActive(false);
	}
}

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

DetectionState ObserverPoint::checkDetection(Candidate *candidate) const {
	double x = candidate->current.getPosition().x;
	if (x > 0) {
		candidate->limitNextStep(x);
		return NOTHING;
	}
	return DETECTED;
}

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

DetectionState ObserverNucleusVeto::checkDetection(Candidate *c) const {
	if (isNucleus(c->current.getId()))
		return NOTHING;
	return VETO;
}

DetectionState ObserverNeutrinoVeto::checkDetection(Candidate *c) const {
	int id = fabs(c->current.getId());
	if ((id == 12) or (id == 14) or (id == 16))
		return NOTHING;
	return VETO;
}

DetectionState ObserverPhotonVeto::checkDetection(Candidate *c) const {
	if (c->current.getId() == 22)
		return NOTHING;
	return VETO;
}

ObserverOutput3D::ObserverOutput3D(std::string fname, bool legacy) :
		legacy(legacy) {
	fout.open(fname.c_str());

	if (legacy) {
		fout << "#CRPropa - Output data file\n";
		fout << "#Format - Particle_Type ";
		fout << "Initial_Particle_Type ";
		fout << "Initial_Position[X,Y,Z](Mpc) ";
		fout << "Initial_Momentum[E(EeV),theta,phi] ";
		fout << "Time(Mpc, light travel distance) ";
		fout << "Position[X,Y,Z](Mpc) ";
		fout << "Momentum[E(EeV),theta,phi]\n";
	} else {
		fout
				<< "# D\tID\tID0\tE\tE0\tX\tY\tZ\tX0\tY0\tZ0\tPx\tPy\tPz\tP0x\tP0y\tP0z\n";
		fout << "#\n";
		fout << "# D           Trajectory length [Mpc]\n";
		fout << "# ID          Particle type (PDG MC numbering scheme)\n";
		fout << "# E           Energy [EeV]\n";
		fout << "# X, Y, Z     Position [Mpc]\n";
		fout << "# Px, Py, Pz  Heading (unit vector of momentum)\n";
		fout << "# Initial state: ID0, E0, ...\n";
		fout << "#\n";
	}
}

ObserverOutput3D::~ObserverOutput3D() {
	fout.close();
}

void ObserverOutput3D::onDetection(Candidate *candidate) const {
	char buffer[256];
	size_t p = 0;

	if (legacy) {
		p += sprintf(buffer + p, "%i ",
				convertToCRPropa2NucleusId(candidate->current.getId()));
		p += sprintf(buffer + p, "%i ",
				convertToCRPropa2NucleusId(candidate->source.getId()));
		Vector3d ipos = candidate->source.getPosition() / Mpc;
		p += sprintf(buffer + p, "%.4f %.4f %.4f ", ipos.x, ipos.y, ipos.z);
		double iPhi = candidate->source.getDirection().getPhi();
		double iTheta = candidate->source.getDirection().getTheta();
		double iE = candidate->source.getEnergy() / EeV;
		p += sprintf(buffer + p, "%.4f %.4f %.4f ", iE, iPhi, iTheta);
		double t = comoving2LightTravelDistance(
				candidate->getTrajectoryLength()) / Mpc;
		p += sprintf(buffer + p, "%.4f ", t);
		Vector3d pos = candidate->current.getPosition() / Mpc;
		p += sprintf(buffer + p, "%.4f %.4f %.4f ", pos.x, pos.y, pos.z);
		double phi = candidate->current.getDirection().getPhi();
		double theta = candidate->current.getDirection().getTheta();
		double E = candidate->current.getEnergy() / EeV;
		p += sprintf(buffer + p, "%.4f %.4f %.4f\n", E, phi, theta);
	} else {
		p += sprintf(buffer + p, "%8.3f\t",
				candidate->getTrajectoryLength() / Mpc);
		p += sprintf(buffer + p, "%10i\t", candidate->current.getId());
		p += sprintf(buffer + p, "%10i\t", candidate->source.getId());
		p += sprintf(buffer + p, "%8.4f\t",
				candidate->current.getEnergy() / EeV);
		p += sprintf(buffer + p, "%8.4f\t",
				candidate->source.getEnergy() / EeV);
		Vector3d pos = candidate->current.getPosition() / Mpc;
		p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", pos.x, pos.y, pos.z);
		Vector3d ipos = candidate->source.getPosition() / Mpc;
		p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", ipos.x, ipos.y,
				ipos.z);
		Vector3d dir = candidate->current.getDirection();
		p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", dir.x, dir.y, dir.z);
		Vector3d idir = candidate->source.getDirection();
		p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\n", idir.x, idir.y,
				idir.z);
	}

#pragma omp critical
	{
		fout.write(buffer, p);
		fout.flush();
	}
}

ObserverOutput1D::ObserverOutput1D(std::string filename, bool legacy) :
		legacy(legacy) {
	fout.open(filename.c_str());

	if (legacy) {
		fout << "#CRPropa - Output data file\n";
		fout << "#Format - Energy(EeV) ";
		fout << "Time(Mpc, light travel distance) ";
		fout << "Initial_Particle_Type ";
		fout << "Initial_Energy(EeV)\n";
	} else {
		fout << "#ID\tE\tD\tID0\tE0\n";
		fout << "#\n";
		fout << "# ID  Particle type\n";
		fout << "# E   Energy [EeV]\n";
		fout << "# D   Comoving trajectory length [Mpc]\n";
		fout << "# ID0 Initial particle type\n";
		fout << "# E0  Initial energy [EeV]\n";
	}
}

ObserverOutput1D::~ObserverOutput1D() {
	fout.close();
}

void ObserverOutput1D::onDetection(Candidate *candidate) const {
	char buffer[256];
	size_t p = 0;

	if (legacy) {
		p += sprintf(buffer + p, "%i ",
				convertToCRPropa2NucleusId(candidate->current.getId()));
		p += sprintf(buffer + p, "%.4f ", candidate->current.getEnergy() / EeV);
		double t = comoving2LightTravelDistance(
				candidate->getTrajectoryLength()) / Mpc;
		p += sprintf(buffer + p, "%.4f ", t);
		p += sprintf(buffer + p, "%i ",
				convertToCRPropa2NucleusId(candidate->source.getId()));
		p += sprintf(buffer + p, "%.4f\n", candidate->source.getEnergy() / EeV);
	} else {
		p += sprintf(buffer + p, "%10i\t", candidate->current.getId());
		p += sprintf(buffer + p, "%8.4f\t",
				candidate->current.getEnergy() / EeV);
		p += sprintf(buffer + p, "%9.4f\t",
				candidate->getTrajectoryLength() / Mpc);
		p += sprintf(buffer + p, "%10i\t", candidate->source.getId());
		p += sprintf(buffer + p, "%8.4f\n",
				candidate->source.getEnergy() / EeV);
	}

#pragma omp critical
	{
		fout.write(buffer, p);
		fout.flush();
	}
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

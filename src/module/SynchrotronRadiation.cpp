#include "crpropa/module/SynchrotronRadiation.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

SynchrotronRadiation::SynchrotronRadiation(ref_ptr<MagneticField> field, bool havePhotons, double limit) {
	Brms = 0.;
	setField(field);
	initSpectrum();
	this->havePhotons = havePhotons;
	this->limit = limit;
	secondaryThreshold = 1e7 * eV;
}

SynchrotronRadiation::SynchrotronRadiation(double Brms, bool havePhotons, double limit) {
	setField(Brms);
	initSpectrum();
	this->havePhotons = havePhotons;
	this->limit = limit;
	secondaryThreshold = 1e7 * eV;
}

void SynchrotronRadiation::setField(ref_ptr<MagneticField> f) {
	this->field = f;
}

void SynchrotronRadiation::setField(double f) {
	this->Brms = f;
}

ref_ptr<MagneticField> SynchrotronRadiation::getField() {
	return field;
}

double SynchrotronRadiation::getBrms() {
	return Brms;
}

void SynchrotronRadiation::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void SynchrotronRadiation::setLimit(double limit) {
	this->limit = limit;
}

void SynchrotronRadiation::setSecondaryThreshold(double t) {
	secondaryThreshold = t;
}

double SynchrotronRadiation::getSecondaryThreshold() const {
	return secondaryThreshold;
}

void SynchrotronRadiation::initSpectrum() {
	std::string filename = getDataPath("synchrotron_spectrum.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"SynchrotronRadiation: could not open file " + filename);

	// clear previously loaded interaction rates
	tabx.clear();
	tabCDF.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabx.push_back(pow(10, a));
				tabCDF.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void SynchrotronRadiation::addPhotons(Candidate *candidate, double loss) const {
	double E = candidate->current.getEnergy();
	double mass = candidate->current.getMass();
	double charge = fabs(candidate->current.getCharge());
	double dE = loss; // energy loss
	double z = candidate->getRedshift();
	double B = 0.;
	if (field.valid()) {
		Vector3d Bvec = field->getField(candidate->current.getPosition());
		B = Bvec.getPerpendicularTo(candidate->current.getDirection()).getR(); // get B field perpendicular to direction of flight
	} else
		B = sqrt(2./3.) * Brms; // represents average vertical component of RMS field strength
	B *= pow(1 + z,2.); // cosmological scaling
	double Ecrit = h_planck / 2. / M_PI * 3./2. * c_light / candidate->current.getMomentum().getR() * charge * B * pow(E / mass / c_squared,3.);

	// draw synchrotron photons as long as their energy is smaller than the energy loss in this propagation step
	Random &random = Random::instance();
	while (dE > 0) {
		// draw random value between 0. and maximum of corresponding cdf
		// choose bin of s where cdf(x) = cdf_rand -> x_rand
		size_t i = random.randBin(tabCDF); // draw random bin (upper bin boundary returned)
		double binWidth = (tabx[i] - tabx[i-1]);
		double x = tabx[i-1] + random.rand() * binWidth; // draw random x uniformly distributed in bin
		double Egamma = x * Ecrit;

		// if the remaining energy is not sufficient check for random accepting
		if (Egamma > dE)
			if (random.rand() > (dE / Egamma))
				break; // not accepted

		// create synchrotron photon and repeat with remaining energy
		dE -= Egamma;
		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(),candidate->current.getPosition());
		if (Egamma > secondaryThreshold) // create only photons with energies above threshold
			candidate->addSecondary(22, Egamma, pos);
	}
}

void SynchrotronRadiation::process(Candidate *candidate) const {
	double charge = fabs(candidate->current.getCharge());
	if (charge == 0)
		return; // only charged particles

	double lf = candidate->current.getLorentzFactor();
	double z = candidate->getRedshift();
	double mass = candidate->current.getMass();
	double gammaBeta = sqrt(pow(lf,2.)-1);
	double E = candidate->current.getEnergy();
	double B = 0.;
	if (field.valid()) {
		Vector3d Bvec = field->getField(candidate->current.getPosition());
		B = Bvec.getPerpendicularTo(candidate->current.getDirection()).getR(); // get B field perpendicular to direction of flight
	} else
		B = sqrt(2./3.) * Brms; // represents average vertical component of RMS field strength
	B *= pow(1 + z,2.); // cosmological scaling
	double dEdx = 2./3./(4. * M_PI * epsilon0) * eplus * eplus * pow(gammaBeta,4.) * pow(charge * B / candidate->current.getMomentum().getR(),2.);
	double Ecrit = h_planck / 2. / M_PI * 3./2. * c_light / candidate->current.getMomentum().getR() * charge * B * pow(E / mass / c_squared,3.);

	double step = candidate->getCurrentStep() / (1 + z); // step size in local frame
	double loss = step * dEdx; // energy loss

	if (havePhotons && Ecrit > secondaryThreshold / 14.) // CDF constant for x > 14 -> no photon energies above threshold possible for Ecrit < Ethr / 14
		addPhotons(candidate, loss);

	if (lf * mass * c_squared - loss <= 0) {
		candidate->setActive(false);
	} else {
		candidate->current.setEnergy(lf * mass * c_squared - loss);
		double losslen = candidate->current.getEnergy() / dEdx;
		candidate->limitNextStep(limit * losslen); // conservative estimate with old dEdx and new E
	}
}

std::string SynchrotronRadiation::getDescription() const {
	std::stringstream s;
	s << "Module for calculation of synchrotron energy loss and creation of synchrotron photons.";
	s << " Have synchrotron photons: " << havePhotons;
	s << ", Energy threshold for production of secondary particles: " << secondaryThreshold;
	if (field.valid())
		s << ", Use Magnetic Field";
	else
		s << ", Use Brms: " << Brms / nG << " nG";
	return s.str();
}

} // namespace crpropa

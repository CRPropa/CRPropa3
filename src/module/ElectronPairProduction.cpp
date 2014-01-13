#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/ParticleID.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

ElectronPairProduction::ElectronPairProduction(PhotonField photonField) :
		photonField(photonField) {
	init();
}

void ElectronPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	init();
}

void ElectronPairProduction::init() {
	switch (photonField) {
	case CMB:
		setDescription("ElectronPairProduction: CMB");
		init(getDataPath("epair_CMB.txt"));
		break;
	case IRB:
		setDescription("ElectronPairProduction: IRB");
		init(getDataPath("epair_IRB.txt"));
		break;
	case CMB_IRB:
		setDescription("ElectronPairProduction: CMB and IRB");
		init(getDataPath("epair_CMB_IRB.txt"));
		break;
	default:
		throw std::runtime_error(
				"ElectronPairProduction: unknown photon background");
	}
}

void ElectronPairProduction::init(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"ElectronPairProduction: could not open file " + filename);

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabLorentzFactor.push_back(a * eV / mass_proton / c_squared);
				tabLossRate.push_back(b / a / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

double ElectronPairProduction::lossLength(int id, double lf, double z) const {
	double Z = chargeNumber(id);
	if (Z == 0)
		return std::numeric_limits<double>::max(); // no pair production on uncharged particles

	lf *= (1 + z);
	if (lf < tabLorentzFactor.front())
		return std::numeric_limits<double>::max(); // below energy threshold

	double rate;
	if (lf < tabLorentzFactor.back())
		rate = interpolate(lf, tabLorentzFactor, tabLossRate); // interpolation
	else
		rate = tabLossRate.back() * pow(lf / tabLorentzFactor.back(), -0.6); // extrapolation

	double A = nucleusMass(id) / mass_proton; // more accurate than massNumber(Id)
	rate *= Z * Z / A * pow(1 + z, 3) * photonFieldScaling(photonField, z);
	return 1. / rate;
}

void ElectronPairProduction::process(Candidate *c) const {
	int id = c->current.getId();
	if (not (isNucleus(id)))
		return; // only nuclei

	double lf = c->current.getLorentzFactor();
	double z = c->getRedshift();
	double step = c->getCurrentStep() / (1 + z); // step size in local frame
	double loss = step / lossLength(id, lf, z);

	c->current.setLorentzFactor(lf * (1 - loss));
}

} // namespace crpropa

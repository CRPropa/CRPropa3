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
				tabEnergy.push_back(a * eV);
				tabLossRate.push_back(b * eV / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

double ElectronPairProduction::lossRate(int id, double E, double z) const {
	double A = massNumberFromNucleusId(id);
	double Z = chargeNumberFromNucleusId(id);

	if (Z < 1)
		return 0; // no pair production on uncharged particles

	double Eeff = E / A * (1 + z);
	if (Eeff < tabEnergy.front())
		return 0; // below energy threshold

	double rate;
	if (Eeff < tabEnergy.back())
		rate = interpolate(Eeff, tabEnergy, tabLossRate);// interpolation
	else
		rate = tabLossRate.back() * pow(Eeff / tabEnergy.back(), 0.4); // extrapolation

	return rate * Z * Z / A * lossRateScaling(photonField, z);
}

void ElectronPairProduction::process(Candidate *c) const {
	int id = c->current.getId();
	if (not(isNucleus(id)))
		return; // this module only handles nucleons and nuclei

	double E = c->current.getEnergy();
	double z = c->getRedshift();
	double step = c->getCurrentStep() / (1 + z); // step size in local frame
	double dE = lossRate(id, E, z) * step;

	c->current.setEnergy(E - dE);
}

} // namespace crpropa

#include "mpc/module/ElectronPairProduction.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace mpc {

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
	// load energy loss rate table
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"ElectronPairProduction: could not open file " + filename);

	std::vector<double> x, y;
	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				energy.push_back(a * eV);
				lossRate.push_back(b * eV / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

void ElectronPairProduction::process(Candidate *candidate) const {
	if (not(candidate->current.isNucleus()))
		return; // this module only handles nucleons / nuclei

	double Z = candidate->current.getChargeNumber();
	double A = candidate->current.getMassNumber();
	double E = candidate->current.getEnergy();
	double z = candidate->getRedshift();
	double EpA = E / A * (1 + z);

	if (Z < 1)
		return; // no electron pair production on uncharged particles

	if (EpA < energy.front())
		return; // below energy threshold

	double rate;
	if (EpA < energy.back())
		rate = interpolate(EpA, energy, lossRate);
	else
		rate = lossRate.back() * pow(EpA / energy.back(), 0.4); // extrapolation

		// step size in local frame
	double step = candidate->getCurrentStep() / (1 + z);

	// dE(E) = Z^2 * loss_rate(E/A) * step
	double dE = Z * Z * rate * photonFieldScaling(photonField, z) * step;
	candidate->current.setEnergy(E - dE);
}

double ElectronPairProduction::energyLossLength(int id, double E) {
	double A = massNumberFromNucleusId(id);
	double Z = chargeNumberFromNucleusId(id);

	if (Z < 1)
		return std::numeric_limits<double>::max();

	double EpA = E / A;
	if (EpA < energy.front())
		return std::numeric_limits<double>::max();

	double rate;
	if (EpA < energy.back())
		rate = interpolate(EpA, energy, lossRate);
	else
		rate = lossRate.back() * pow(EpA / energy.back(), 0.4);

	double lossRate = Z * Z * rate / E;
	return 1. / lossRate;
}

} // namespace mpc

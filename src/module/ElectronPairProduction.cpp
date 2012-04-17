#include "mpc/module/ElectronPairProduction.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace mpc {

ElectronPairProduction::ElectronPairProduction(int photonField) {
	init(photonField);
}

void ElectronPairProduction::init(int photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("ElectronPairProduction:CMB");
		init(getDataPath("ElectronPairProduction/cmb.txt"));
		break;
	case IRB:
		setDescription("ElectronPairProduction:IRB");
		init(getDataPath("ElectronPairProduction/ir.txt"));
		break;
	case CMB_IRB:
		setDescription("ElectronPairProduction:CMB_IRB");
		init(getDataPath("ElectronPairProduction/cmbir.txt"));
		break;
	default:
		throw std::runtime_error("mpc::ElectronPairProduction: unknown photon background");
	}
}

void ElectronPairProduction::init(std::string filename) {
	// load energy loss rate table
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"mpc::ElectronPairProduction: could not open file " + filename);

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			x.push_back(a * eV);
			y.push_back(b * eV / Mpc);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

void ElectronPairProduction::process(Candidate *candidate) const {
	double A = candidate->current.getMassNumber();
	double E = candidate->current.getEnergy();
	double z = candidate->getRedshift();

	// energy per nucleon
	double xi = E / A * (1 + z);

	// no data for EpA < 10^15 eV
	if (xi < 1.e15 * eV)
		return;

	double rate;
	if (xi < 1.e22 * eV) {
		// index of next lower energy bin
		int i = floor(log10(xi / x[0]) * 69 / 7);
		// linear interpolation: y(x) = y0 + dy/dx * (x-x0)
		rate = y[i] + (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (xi - x[i]);
	} else {
		// extrapolation for higher energies
		rate = y[69] * pow(xi / x[69], 0.4);
	}

	// dE(E) = Z^2 * loss_rate(E/A) * step
	double step = candidate->getCurrentStep();
	double Z = candidate->current.getChargeNumber();
	double dE = Z * Z * rate * pow(1 + z, 3) * step;
	candidate->current.setEnergy(E - dE);
}

} // namespace mpc

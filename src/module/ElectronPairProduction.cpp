#include "mpc/module/ElectronPairProduction.h"

#include <fstream>
#include <vector>
#include <limits>
#include <stdexcept>

namespace mpc {

ElectronPairProduction::ElectronPairProduction(int photonField) {
	init(photonField);
}

ElectronPairProduction::ElectronPairProduction() {
	init(CMBIR);
}

void ElectronPairProduction::init(int photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		init(getDataPath("ElectronPairProduction/cmb.txt"));
		break;
	case IR:
		init(getDataPath("ElectronPairProduction/ir.txt"));
		break;
	case CMBIR:
		init(getDataPath("ElectronPairProduction/cmbir.txt"));
		break;
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
			if (infile) {
				x.push_back(a * eV);
				y.push_back(b * eV / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

std::string ElectronPairProduction::getDescription() const {
	switch (photonField) {
	case CMB: {
		return "Electron-pair production on CMB";
		break;
	}
	case IR: {
		return "Electron-pair production on IB";
		break;
	}
	case CMBIR: {
		return "Electron-pair production on CMB + IR";
		break;
	}
	default: {
		return "Electron-pair production (unknown)";
		break;
	}
	}
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

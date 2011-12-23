#include "mpc/module/ElectronPairProduction.h"

#include <fstream>
#include <vector>
#include <limits>

namespace mpc {

ElectronPairProduction::ElectronPairProduction(PhotonField photonField) {
	init(photonField);
}

ElectronPairProduction::ElectronPairProduction() {
	init(CMBIR);
}

void ElectronPairProduction::init(PhotonField photonField) {
	field = photonField;
	switch (photonField) {
	case CMB:
		init("data/ElectronPairProduction/cmb.txt");
		break;
	case IR:
		init("data/ElectronPairProduction/ir.txt");
		break;
	case CMBIR:
		init("data/ElectronPairProduction/cmbir.txt");
		break;
	}
}

void ElectronPairProduction::init(std::string filename) {
	// load energy loss rate table
	std::ifstream infile(filename.c_str());
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
	static std::string u("Electron-pair production (unknown)");

	switch (field) {
	case CMB: {
		static std::string s("Electron-pair production (CMB)");
		return s;
		break;
	}
	case IR: {
		static std::string s("Electron-pair production (CMB)");
		return s;
		break;
	}
	case CMBIR: {
		static std::string s("Electron-pair production (CMB)");
		return s;
		break;
	}
	default: {
		return u;
		break;
	}
	}
	return u;
}

void ElectronPairProduction::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	double A = candidate->current.getMassNumber();
	double E = candidate->current.getEnergy();
	double z = candidate->getRedshift();

	// energy per nucleon for lookup in proton table
	double EpA = E / A * (1 + z);

	// no data for EpA < 10^15 eV
	if (EpA < 1.e15 * eV)
		return;

	double rate;
	if (EpA < 1.e22 * eV) {
		// index of next lower energy bin
		int i = floor(log10(EpA / x[0]) * 69 / 7);
		// linear interpolation: y(x) = y0 + dy/dx * (x-x0)
		rate = y[i] + (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (EpA - x[i]);
	} else {
		// extrapolation for higher energies
		rate = y[69] * pow(EpA / x[69], 0.4);
	}

	// dE(E) = Z^2 * loss_rate(E/A) * step
	double step = candidate->getCurrentStep();
	double Z = candidate->current.getChargeNumber();
	double dE = Z * Z * rate * pow(1 + z, 3) * step;
	candidate->current.setEnergy(E - dE);
}

} // namespace mpc

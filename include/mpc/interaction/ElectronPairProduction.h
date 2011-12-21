#ifndef ELECTRONPAIRPRODUCTION_H_
#define ELECTRONPAIRPRODUCTION_H_

#include "mpc/Module.h"
#include "mpc/ParticleState.h"
#include <fstream>
#include <vector>

namespace mpc {

class ElectronPairProduction: public Module {
private:
	std::vector<double> y; // energy loss rate table for protons in [J/m]
	std::vector<double> x; // energy table in [J]

public:
	enum PhotonField {
		CMB, IR, CMBIR
	};

	ElectronPairProduction(PhotonField photonField) {
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

	ElectronPairProduction() {
		init("data/ElectronPairProduction/cmbir.txt");
	}

	void init(std::string filename) {
		// load energy loss rate table
		std::ifstream infile(filename.c_str());
		char header[256];
		for (int i = 0; i < 3; i++)
			infile.getline(header, 255);
		double a, b;
		while (infile.good()) {
			infile >> a >> b;
			x.push_back(a * eV);
			y.push_back(b * eV / Mpc);
		}
		infile.close();
	}

	std::string getDescription() const {
		return "Electron-pair production";
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		double E = candidate->current.getEnergy();
		double A = candidate->current.getMassNumber();
		double z = candidate->getRedshift();

		// redshift modified energy per nucleon for lookup in proton table
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
};

} // namespace mpc

#endif /* ELECTRONPAIRPRODUCTION_H_ */

#ifndef ELECTRONPAIRPRODUCTION_H_
#define ELECTRONPAIRPRODUCTION_H_

#include "mpc/Module.h"
#include "mpc/ParticleState.h"
#include <fstream>

namespace mpc {

class ElectronPairProduction: public Module {
private:
	double y[69]; // tabulated electron pair production energy loss rate for protons in [J / m]
	double x[69]; // tabulated energy in [J]
	double dlx; // value of logarithmic spacing

public:
	ElectronPairProduction() {
		// load energy loss rate table
		std::ifstream infile("data/pair_rate_cmbir.table");
		char str[256];
		for (int i = 0; i < 3; i++)
			infile.getline(str, 255); // skip header
		double a, b;
		for (int i = 0; i < 70; i++) {
			infile >> a >> b;
			x[i] = a * eV; // convert [eV] -> [J]
			y[i] = b * eV / Mpc; // convert [eV / Mpc] -> [J / m]
		}
		infile.close();
		dlx = log10(x[69] / x[0]) / 69;
	}

	std::string getDescription() const {
		return "Electron-pair production";
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		double step = candidate->getCurrentStep();
		double E = candidate->current.getEnergy();
		double Z = candidate->current.getChargeNumber();
		double mass = candidate->current.getMass();

		// index of next lower energy bin
		int i = floor(log10(E / x[0]) / dlx);

		// energy loss for nucleus of mass M
		// dE(E) = m_p / m * Z^2 * loss_rate(E) * step
		// linear interpolation: y(x) = y0 + dy/dx * (x-x0)
		double rate = y[i] + (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (E - x[i]);
		double dE = Z * Z * mass_proton / mass * rate * step;

		std::cout << E / EeV << std::endl;
		std::cout << dE / EeV << std::endl << std::endl;

		// update particle
		candidate->current.setEnergy(E - dE);
	}
};

} // namespace mpc

#endif /* ELECTRONPAIRPRODUCTION_H_ */

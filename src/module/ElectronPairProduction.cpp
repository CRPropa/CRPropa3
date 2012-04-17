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
		throw std::runtime_error(
				"mpc::ElectronPairProduction: unknown photon background");
	}
}

void ElectronPairProduction::init(std::string filename) {
	// load energy loss rate table
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error(
				"mpc::ElectronPairProduction: could not open file " + filename);

	std::vector<double> x, y;
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

	acc = gsl_interp_accel_alloc();
	lossRate = gsl_spline_alloc(gsl_interp_linear, x.size());
	gsl_spline_init(lossRate, &x[0], &y[0], x.size());
	xMin = x.front();
	xMax = x.back();
	yMax = y.back();
}

void ElectronPairProduction::process(Candidate *candidate) const {
	double A = candidate->current.getMassNumber();
	double E = candidate->current.getEnergy();
	double z = candidate->getRedshift();

	double EpA = E / A * (1 + z);
	if (EpA < xMin)
		return;

	double rate;
	if (EpA < xMax)
		rate = gsl_spline_eval(lossRate, EpA, acc);
	else
		// extrapolation for higher energies
		rate = yMax * pow(EpA / xMax, 0.4);

	// dE(E) = Z^2 * loss_rate(E/A) * step
	double step = candidate->getCurrentStep();
	double Z = candidate->current.getChargeNumber();
	double dE = Z * Z * rate * pow(1 + z, 3) * step;
	candidate->current.setEnergy(E - dE);
}

} // namespace mpc

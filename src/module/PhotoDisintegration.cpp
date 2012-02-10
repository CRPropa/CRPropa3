#include "mpc/module/PhotoDisintegration.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>

namespace mpc {

PhotoDisintegration::PhotoDisintegration() {
	// create Lorentz factor sample points
	std::vector<double> x;
	for (size_t i = 0; i < 200; i++) {
		x.push_back(6.0 + i * 8.0 / 199);
	}

	// load disintegration data
	std::ifstream infile("data/PhotoDisintegration/rate.txt");
	std::string line, cell;
	while (std::getline(infile, line)) {
		std::stringstream lineStream(line);

		std::getline(lineStream, cell, ' ');
		int id = atoi(cell.c_str()); // nucleus id

		DisintegrationMode mode;
		std::getline(lineStream, cell, ' ');
		mode.channel = atoi(cell.c_str()); // disintegration channel

		std::vector<double> y;
		while (std::getline(lineStream, cell, ' ')) {
			double a = atof(cell.c_str()); // disintegration rate
			y.push_back(a / Mpc);
		}
		mode.rate = gsl_spline_alloc(gsl_interp_linear, 200);
		gsl_spline_init(mode.rate, &x[0], &y[0], 200);

		modeMap[id].push_back(mode);
	}
	infile.close();
	acc = gsl_interp_accel_alloc();
}

std::string PhotoDisintegration::getDescription() const {
	return "Photo-Disintegration";
}

void PhotoDisintegration::process(Candidate *candidate) {
	double step = candidate->getCurrentStep();

	while (true) {
		int id = candidate->current.getId();

		// set a new interaction if necessary
		if (id != cached_id) {
			// return if no disintegration data
			if (setNextInteraction(candidate) == false)
				return;
			cached_id = id;
		}
		// if counter not over, reduce and return
		if (cached_distance > step) {
			cached_distance -= step;
			candidate->limitNextStep(cached_distance);
			return;
		}
		// counter over: disintegrate
		cached_id = 0;
		step -= cached_distance;
		performInteraction(candidate);
	}
}

bool PhotoDisintegration::setNextInteraction(Candidate *candidate) {
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double lg = log10(candidate->current.getLorentzFactor() * (1 + z));

	// no disintegration data
	if (modeMap[id].size() == 0)
		return false;

	// out of energy range
	if ((lg < 6) or (lg > 14))
		return false;

	// find channel with minimum random decay distance
	cached_distance = std::numeric_limits<double>::max();
	for (size_t i = 0; i < modeMap[id].size(); i++) {
		double rate = gsl_spline_eval(modeMap[id][i].rate, lg, acc);
		double d = -log(rand.rand()) / rate;
		if (d > cached_distance)
			continue;
		cached_distance = d;
		cached_channel = modeMap[id][i].channel;
	}

	cached_distance /= pow((1 + z), 3); // CMB density increases with (1+z)^3
	return true;
}

void PhotoDisintegration::performInteraction(Candidate *candidate) {
	// parse disintegration channel
	int nNeutron = cached_channel / 100000;
	int nProton = (cached_channel % 100000) / 10000;
	int nH2 = (cached_channel % 10000) / 1000;
	int nH3 = (cached_channel % 1000) / 100;
	int nHe3 = (cached_channel % 100) / 10;
	int nHe4 = (cached_channel % 10);

	int dA = -nNeutron - nProton - 2 * nH2 - 3 * nH3 - 3 * nHe3 - 4 * nHe4;
	int dZ = -nProton - nH2 - nH3 - 2 * nHe3 - 2 * nHe4;

	int id = candidate->current.getId();
	int A = getMassNumberFromNucleusId(id);
	int Z = getChargeNumberFromNucleusId(id);
	double EpA = candidate->current.getEnergy() / double(A);

	// update particle
	candidate->current.setId(getNucleusId(A + dA, Z + dZ));
	candidate->current.setEnergy(EpA * (A + dA));

	// create secondaries
	for (size_t i = 0; i < nNeutron; i++) {
		candidate->addSecondary(getNucleusId(1, 0), EpA);
	}
	for (size_t i = 0; i < nProton; i++) {
		candidate->addSecondary(getNucleusId(1, 1), EpA);
	}
	for (size_t i = 0; i < nH2; i++) {
		candidate->addSecondary(getNucleusId(2, 1), EpA * 2);
	}
	for (size_t i = 0; i < nH3; i++) {
		candidate->addSecondary(getNucleusId(3, 1), EpA * 3);
	}
	for (size_t i = 0; i < nHe3; i++) {
		candidate->addSecondary(getNucleusId(3, 2), EpA * 3);
	}
	for (size_t i = 0; i < nHe4; i++) {
		candidate->addSecondary(getNucleusId(4, 2), EpA * 4);
	}
}

} // namespace mpc

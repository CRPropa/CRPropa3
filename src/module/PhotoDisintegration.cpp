#include "mpc/module/PhotoDisintegration.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>

namespace mpc {

PhotoDisintegration::PhotoDisintegration() {
	cached_id = 0;
	// load disintegration data
	std::ifstream infile("data/PhotoDisintegration/meanFreePath.txt");
	std::string line, cell;
	while (std::getline(infile, line)) {
		std::stringstream lineStream(line);
		std::getline(lineStream, cell, ' ');
		int id = atoi(cell.c_str()); // nucleus id
		std::getline(lineStream, cell, ' ');
		DisintegrationMode mode;
		mode.channel = atoi(cell.c_str()); // disintegration channel
		while (std::getline(lineStream, cell, ' ')) {
			double mfp = 1 / atof(cell.c_str()) / Mpc; // mean free path [m]
			mode.y.push_back(mfp);
		}
		modeMap[id].push_back(mode);
	}
	infile.close();
	for (size_t i = 0; i < 200; i++) {
		x.push_back(6.0 + i * 8.0 / 199); // log(gamma)
	}
}

std::string PhotoDisintegration::getDescription() const {
	return "Photo-Disintegration";
}

void PhotoDisintegration::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	int id = candidate->current.getId();
	std::vector<DisintegrationMode> modes = modeMap[id];

	// no disintegration modes for this nucleus
	if (modes.size() == 0)
		return;

	double gamma = candidate->current.getLorentzFactor();
	double z = candidate->getRedshift();
	double lg = log10(gamma * (1 + z));

	// no data for log10(gamma) < 6 or log10(gamma) > 14
	if (lg < 6)
		return;
	if (lg > 14)
		return;

	// check if a disintegration is already cached, if not set it
	if (id != cached_id) {
		cached_id = id;
		// find decay mode with minimum random decay distance
		cached_distance = std::numeric_limits<double>::max();
		for (size_t i = 0; i < modes.size(); i++) {
			double mfp = getMeanFreePath(modes[i].y, lg);
			double d = -log(mtrand.rand()) * mfp;
			if (d > cached_distance)
				continue;
			cached_distance = d;
			cached_channel = modes[i].channel;
		}
	}

	// check if life-time is over, reduce life-time
	double step = candidate->getCurrentStep();
	if (cached_distance > step) {
		cached_distance -= step;
		candidate->limitNextStep(cached_distance * gamma);
		return;
	}

	// else life time is over, decay
	cached_id = 0;
	disintegrate(candidate, secondaries);
}

double PhotoDisintegration::getMeanFreePath(std::vector<double> y,
		double lg) {
	// index of next lower gamma bin
	size_t i = (lg - 6) / 200;
	// linear interpolation: y(x) = y0 + dy/dx * (x-x0)
	return y[i] + (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (lg - x[i]);
}

void PhotoDisintegration::disintegrate(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
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

	candidate->current.setId(getNucleusId(A + dA, Z + dZ));
	candidate->current.setEnergy(EpA * (A + dA));

	for (size_t i = 0; i < nNeutron; i++) {
		createSecondary(candidate, secondaries, getNucleusId(1, 0), EpA);
	}
	for (size_t i = 0; i < nProton; i++) {
		createSecondary(candidate, secondaries, getNucleusId(1, 1), EpA);
	}
	for (size_t i = 0; i < nH2; i++) {
		createSecondary(candidate, secondaries, getNucleusId(2, 1), EpA * 2);
	}
	for (size_t i = 0; i < nH3; i++) {
		createSecondary(candidate, secondaries, getNucleusId(3, 1), EpA * 3);
	}
	for (size_t i = 0; i < nHe3; i++) {
		createSecondary(candidate, secondaries, getNucleusId(3, 2), EpA * 3);
	}
	for (size_t i = 0; i < nHe4; i++) {
		createSecondary(candidate, secondaries, getNucleusId(4, 2), EpA * 4);
	}
}

void PhotoDisintegration::createSecondary(Candidate *candidate,
		std::vector<Candidate *> &secondaries, int id, double energy) {
	ParticleState initial = candidate->current;
	initial.setEnergy(energy);
	initial.setId(id);
	Candidate secondary;
	secondary.current = initial;
	secondary.initial = initial;
	secondary.setNextStep(candidate->getCurrentStep());
	secondaries.push_back(&secondary);
}

} // namespace mpc

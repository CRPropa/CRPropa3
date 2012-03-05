#include "mpc/module/PhotoDisintegration.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>

namespace mpc {

PhotoDisintegration::PhotoDisintegration() {
	name = "mpc::PhotoDisintegration";

	// create spline x-axis
	std::vector<double> x;
	for (size_t i = 0; i < 200; i++) {
		x.push_back(6.0 + i * 8.0 / 199);
	}

	// load disintegration table
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

		PDTable[id].push_back(mode);
	}
	infile.close();
	acc = gsl_interp_accel_alloc();
}

std::string PhotoDisintegration::getDescription() const {
	return "Photo-disintegration";
}

void PhotoDisintegration::process(Candidate *candidate) {
	double step = candidate->getCurrentStep();
	InteractionState interaction;

	while (true) {
		// set a new interaction if necessary
		if (not (candidate->getInteractionState(name, interaction))) {
			if (not (setNextInteraction(candidate)))
				return;
		}

		// get the interaction state
		candidate->getInteractionState(name, interaction);

		// if not over, reduce distance and return
		if (interaction.distance > step) {
			interaction.distance -= step;
			candidate->limitNextStep(interaction.distance);
			candidate->setInteractionState(name, interaction);
			return;
		}

		// counter over: interact
		step -= interaction.distance;
		performInteraction(candidate);
	}
}

bool PhotoDisintegration::setNextInteraction(Candidate *candidate) {
	int id = candidate->current.getId();

	// check if disintegration data available
	if (PDTable[id].size() == 0)
		return false;

	double z = candidate->getRedshift();
	double lg = log10(candidate->current.getLorentzFactor() * (1 + z));

	// check if out of energy range
	if ((lg < 6) or (lg > 14))
		return false;

	// find channel with minimum random decay distance
	InteractionState interaction;
	interaction.distance = std::numeric_limits<double>::max();
	for (size_t i = 0; i < PDTable[id].size(); i++) {
		double rate = gsl_spline_eval(PDTable[id][i].rate, lg, acc);
		double d = -log(random.rand()) / rate;
		if (d > interaction.distance)
			continue;
		interaction.distance = d;
		interaction.channel = PDTable[id][i].channel;
	}

	// CMB density increases with (1+z)^3 -> free distance decreases accordingly
	interaction.distance /= pow((1 + z), 3);

	candidate->setInteractionState(name, interaction);
	return true;
}

void PhotoDisintegration::performInteraction(Candidate *candidate) {
	InteractionState interaction;
	candidate->getInteractionState(name, interaction);
	candidate->clearInteractionStates();

	// parse disintegration channel
	int nNeutron = interaction.channel / 100000;
	int nProton = (interaction.channel % 100000) / 10000;
	int nH2 = (interaction.channel % 10000) / 1000;
	int nH3 = (interaction.channel % 1000) / 100;
	int nHe3 = (interaction.channel % 100) / 10;
	int nHe4 = (interaction.channel % 10);

	int dA = -nNeutron - nProton - 2 * nH2 - 3 * nH3 - 3 * nHe3 - 4 * nHe4;
	int dZ = -nProton - nH2 - nH3 - 2 * nHe3 - 2 * nHe4;

	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	double EpA = candidate->current.getEnergy() / double(A);

	// update particle
	candidate->current.setId(getNucleusId(A + dA, Z + dZ));
	candidate->current.setEnergy(EpA * (A + dA));

	// create secondaries
	for (size_t i = 0; i < nNeutron; i++)
		candidate->addSecondary(getNucleusId(1, 0), EpA);
	for (size_t i = 0; i < nProton; i++)
		candidate->addSecondary(getNucleusId(1, 1), EpA);
	for (size_t i = 0; i < nH2; i++)
		candidate->addSecondary(getNucleusId(2, 1), EpA * 2);
	for (size_t i = 0; i < nH3; i++)
		candidate->addSecondary(getNucleusId(3, 1), EpA * 3);
	for (size_t i = 0; i < nHe3; i++)
		candidate->addSecondary(getNucleusId(3, 2), EpA * 3);
	for (size_t i = 0; i < nHe4; i++)
		candidate->addSecondary(getNucleusId(4, 2), EpA * 4);
}

} // namespace mpc

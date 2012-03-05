#include "mpc/module/PhotoDisintegration.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdexcept>

#include <kiss/convert.h>

namespace mpc {

enum _Constants {
	SAMPLE_COUNT = 200
};

PhotoDisintegration::PhotoDisintegration() {
	name = "mpc::PhotoDisintegration";
	acc = gsl_interp_accel_alloc();

	// create spline x-axis
	std::vector<double> x(SAMPLE_COUNT);
	for (size_t i = 0; i < 200; i++) {
		x[i] = 6.0 + i * 8.0 / (SAMPLE_COUNT - 1);
	}

	std::vector<double> y(SAMPLE_COUNT);

	// load disintegration table
	std::ifstream infile("data/PhotoDisintegration/rate.txt");
	std::string line, cell;
	size_t lineNo = 0;
	while (std::getline(infile, line)) {
		lineNo++;
		std::stringstream lineStream(line);

		int id = 0;
		lineStream >> id;

		DisintegrationMode mode;
		lineStream >> mode.channel; // disintegration channel

		for (size_t i = 0; i < SAMPLE_COUNT; i++) {
			lineStream >> y[i];
			y[i] /= Mpc;

			if (!lineStream)
				throw std::runtime_error(
						"mpc::PhotoDisintegration: not enough entries in rate.txt. Expected 200 got "
								+ kiss::str(i) + " in line "
								+ kiss::str(lineNo));
		}

		mode.rate = gsl_spline_alloc(gsl_interp_linear, SAMPLE_COUNT);
		gsl_spline_init(mode.rate, &x[0], &y[0], SAMPLE_COUNT);

		PDTable[id].push_back(mode);
	}
}

PhotoDisintegration::~PhotoDisintegration() {
	gsl_interp_accel_free(acc);
	DisintegrationModeMap::iterator iModes;
	for (iModes = PDTable.begin(); iModes != PDTable.end(); iModes++) {
		for (size_t iMode = 0; iMode < iModes->second.size(); iMode++)
			gsl_spline_free(iModes->second[iMode].rate);
	}
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

	DisintegrationModeMap::iterator iMode = PDTable.find(id);
	if (iMode == PDTable.end())
		return false;

	std::vector<DisintegrationMode> &modes = iMode->second;

	// check if disintegration data available
	if (modes.size() == 0)
		return false;

	double z = candidate->getRedshift();
	double lg = log10(candidate->current.getLorentzFactor() * (1 + z));

	// check if out of energy range
	if ((lg < 6) or (lg > 14))
		return false;

	// find channel with minimum random decay distance
	InteractionState interaction;
	interaction.distance = std::numeric_limits<double>::max();
	for (size_t i = 0; i < modes.size(); i++) {
		double rate = gsl_spline_eval(modes[i].rate, lg, acc);
		double d = -log(random.rand()) / rate;
		if (d > interaction.distance)
			continue;
		interaction.distance = d;
		interaction.channel = modes[i].channel;
	}

	// CMB density increases with (1+z)^3 -> free distance decreases accordingly
	interaction.distance /= pow((1 + z), 3);

	candidate->setInteractionState(name, interaction);
	return true;
}

inline int digit(const int &value, const int &d) {
	return (value % d * 10) / d;
}

void PhotoDisintegration::performInteraction(Candidate *candidate) {
	InteractionState interaction;
	candidate->getInteractionState(name, interaction);
	candidate->clearInteractionStates();

	// parse disintegration channel
	int nNeutron = digit(interaction.channel, 100000);
	int nProton = digit(interaction.channel, 10000);
	int nH2 = digit(interaction.channel, 1000);
	int nH3 = digit(interaction.channel, 100);
	int nHe3 = digit(interaction.channel, 10);
	int nHe4 = digit(interaction.channel, 1);

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

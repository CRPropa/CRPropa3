#include "mpc/module/PhotoDisintegration.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include <kiss/convert.h>

namespace mpc {

enum _Constants {
	SAMPLE_COUNT = 200
};

PhotoDisintegration::PhotoDisintegration(int photonField) {
	init(photonField);
}

void PhotoDisintegration::init(int photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("PhotoDisintegration:CMB");
		init(getDataPath("PhotoDisintegration/PDtable_CMB.txt"));
		break;
	case IRB:
		setDescription("PhotoDisintegration:IRB");
		init(getDataPath("PhotoDisintegration/PDtable_IRB.txt"));
		break;
	case CMB_IRB:
		setDescription("PhotoDisintegration:CMB_IRB");
		init(getDataPath("PhotoDisintegration/PDtable_CMB_IRB.txt"));
		break;
	default:
		throw std::runtime_error(
				"mpc::PhotoDisintegration: unknown photon background");
	}
}

void PhotoDisintegration::init(std::string filename) {
	std::ifstream infile(filename.c_str());

	acc = gsl_interp_accel_alloc();

	// create spline x-axis
	std::vector<double> x(SAMPLE_COUNT);
	for (size_t i = 0; i < 200; i++)
		x[i] = 6.0 + i * 8.0 / (SAMPLE_COUNT - 1);

	std::vector<double> y(SAMPLE_COUNT);

	if (!infile.good())
		throw std::runtime_error(
				"mpc::PhotoDisintegration: could not open file " + filename);

	std::string line;
	size_t lineNo = 0;
	while (std::getline(infile, line)) {
		lineNo++;
		if (line[0] == '#')
			continue;
		std::stringstream lineStream(line);

		int Z, A;
		lineStream >> Z; // charge number
		lineStream >> A; // mass number

		DisintegrationMode mode;
		lineStream >> mode.channel; // disintegration channel

		for (size_t i = 0; i < SAMPLE_COUNT; i++) {
			lineStream >> y[i];
			y[i] /= Mpc; // disintegration rate in [1/m]

			if (!lineStream)
				throw std::runtime_error(
						"mpc::PhotoDisintegration: Not enough entries in rate.txt. Expected 200 got "
								+ kiss::str(i) + " in line "
								+ kiss::str(lineNo));
		}

		mode.rate = gsl_spline_alloc(gsl_interp_linear, SAMPLE_COUNT);
		gsl_spline_init(mode.rate, &x[0], &y[0], SAMPLE_COUNT);

		PDTable[getNucleusId(A, Z)].push_back(mode);
	}

	infile.close();
}

PhotoDisintegration::~PhotoDisintegration() {
	gsl_interp_accel_free(acc);
	DisintegrationModeMap::iterator iModes;
	for (iModes = PDTable.begin(); iModes != PDTable.end(); iModes++) {
		for (size_t iMode = 0; iMode < iModes->second.size(); iMode++)
			gsl_spline_free(iModes->second[iMode].rate);
	}
}

bool PhotoDisintegration::setNextInteraction(Candidate *candidate,
		InteractionState &interaction) const {
	int id = candidate->current.getId();

	// check if disintegration data available
	DisintegrationModeMap::const_iterator iMode = PDTable.find(id);
	if (iMode == PDTable.end())
		return false;

	const std::vector<DisintegrationMode> &modes = iMode->second;

	// CMB energy increases with (1+z), increase nucleus energy accordingly
	double z = candidate->getRedshift();
	double lg = log10(candidate->current.getLorentzFactor() * (1 + z));

	// check if out of energy range
	if ((lg < 6) or (lg > 14))
		return false;

	// find channel with minimum random decay distance
	interaction.distance = std::numeric_limits<double>::max();
	for (size_t i = 0; i < modes.size(); i++) {
		double rate = gsl_spline_eval(modes[i].rate, lg, acc);
		double d = -log(Random::instance().rand()) / rate;
		if (d > interaction.distance)
			continue;
		interaction.distance = d;
		interaction.channel = modes[i].channel;
	}

	// CMB density increases with (1+z)^3 -> free distance decreases accordingly
	interaction.distance /= pow((1 + z), 3);

	candidate->setInteractionState(getDescription(), interaction);
	return true;
}

void PhotoDisintegration::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(getDescription(), interaction);
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
	int nA = A + dA;
	if (nA > 0) {
		candidate->current.setId(getNucleusId(A + dA, Z + dZ));
		candidate->current.setEnergy(EpA * (A + dA));
	} else {
		candidate->setActive(false);
	}

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

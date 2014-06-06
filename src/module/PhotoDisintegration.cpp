#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {

PhotoDisintegration::PhotoDisintegration(PhotonField photonField, double l) {
	limit = l;
	init(photonField);
}

void PhotoDisintegration::setLimit(double l) {
	limit = l;
}

void PhotoDisintegration::init(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("PhotoDisintegration: CMB");
		init(getDataPath("photodis_CMB.txt"));
		break;
	case IRB:
		setDescription("PhotoDisintegration: IRB");
		init(getDataPath("photodis_IRB.txt"));
		break;
	default:
		throw std::runtime_error(
				"crpropa::PhotoDisintegration: unknown photon background");
	}
}

void PhotoDisintegration::init(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error(
				"crpropa::PhotoDisintegration: could not open file "
						+ filename);

	lgmin = 6;
	lgmax = 13.96;

	// clear previously loaded interaction rates
	pdTable.clear();
	pdTable.resize(27 * 31);

	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;
		std::stringstream lineStream(line);

		int Z, N;
		lineStream >> Z;
		lineStream >> N;

		PDMode pd;
		lineStream >> pd.channel;

		double r = 0;
		for (size_t i = 0; i < 200; i++) {
			lineStream >> r;
			pd.rate.push_back(r / Mpc);
		}

		pdTable[Z * 31 + N].push_back(pd);
	}

	infile.close();
}

void PhotoDisintegration::process(Candidate *candidate) const {
	// the loop will be executed at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (not (isNucleus(id)))
			return;

		int A = massNumber(id);
		int Z = chargeNumber(id);
		int N = A - Z;

		// check if disintegration data available
		if ((Z > 26) or (N > 30))
			return;
		std::vector<PDMode> pdModes = pdTable[Z * 31 + N];
		if (pdModes.size() == 0)
			return;

		// check if in tabulated energy range
		double z = candidate->getRedshift();
		double lg = log10(candidate->current.getLorentzFactor() * (1 + z));
		if ((lg <= lgmin) or (lg >= lgmax))
			return;

		// find disintegration channel with minimum random decay distance
		Random &random = Random::instance();
		double randDistance = std::numeric_limits<double>::max();
		int channel;
		double totalRate = 0;

		// comological scaling of interaction distance (comoving units)
		double scaling = pow(1 + z, 2) * photonFieldScaling(photonField, z);

		for (size_t i = 0; i < pdModes.size(); i++) {
			double rate = interpolateEquidistant(lg, lgmin, lgmax,
					pdModes[i].rate);
			rate *= scaling;
			totalRate += rate;
			double d = -log(random.rand()) / rate;
			if (d > randDistance)
				continue;
			randDistance = d;
			channel = pdModes[i].channel;
		}

		// check if interaction doesn't happen
		if (step < randDistance) {
			// limit next step to a fraction of the mean free path
			candidate->limitNextStep(limit / totalRate);
			return;
		}

		// interact and repeat with remaining step
		performInteraction(candidate, channel);
		step -= randDistance;
	} while (step > 0);
}

void PhotoDisintegration::performInteraction(Candidate *candidate,
		int channel) const {
	// interpret disintegration channel
	int nNeutron = digit(channel, 100000);
	int nProton = digit(channel, 10000);
	int nH2 = digit(channel, 1000);
	int nH3 = digit(channel, 100);
	int nHe3 = digit(channel, 10);
	int nHe4 = digit(channel, 1);

	int dA = -nNeutron - nProton - 2 * nH2 - 3 * nH3 - 3 * nHe3 - 4 * nHe4;
	int dZ = -nProton - nH2 - nH3 - 2 * nHe3 - 2 * nHe4;

	int id = candidate->current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);
	double EpA = candidate->current.getEnergy() / A;

	// update particle
	int nA = A + dA;
	if (nA > 0) {
		candidate->current.setId(nucleusId(A + dA, Z + dZ));
		candidate->current.setEnergy(EpA * (A + dA));
	} else {
		candidate->setActive(false);
	}

	// create secondaries
	for (size_t i = 0; i < nNeutron; i++)
		candidate->addSecondary(nucleusId(1, 0), EpA);
	for (size_t i = 0; i < nProton; i++)
		candidate->addSecondary(nucleusId(1, 1), EpA);
	for (size_t i = 0; i < nH2; i++)
		candidate->addSecondary(nucleusId(2, 1), EpA * 2);
	for (size_t i = 0; i < nH3; i++)
		candidate->addSecondary(nucleusId(3, 1), EpA * 3);
	for (size_t i = 0; i < nHe3; i++)
		candidate->addSecondary(nucleusId(3, 2), EpA * 3);
	for (size_t i = 0; i < nHe4; i++)
		candidate->addSecondary(nucleusId(4, 2), EpA * 4);
}

double PhotoDisintegration::lossLength(int id, double E, double z) {
	// check if nucleus
	if (not (isNucleus(id)))
		return 0;

	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;

	// check if disintegration data available
	if ((Z > 26) or (N > 30))
		return std::numeric_limits<double>::max();
	std::vector<PDMode> pdModes = pdTable[Z * 31 + N];
	if (pdModes.size() == 0)
		return std::numeric_limits<double>::max();

	// check if in tabulated energy range
	double lg = log10(E / (nucleusMass(id) * c_squared)) * (1 + z);
	if ((lg <= lgmin) or (lg >= lgmax))
		return std::numeric_limits<double>::max();

	// total rate from all disintegration channels
	double lossRate = 0;
	for (size_t i = 0; i < pdModes.size(); i++) {
		int dA = 0;
		dA += 1 * digit(pdModes[i].channel, 100000);
		dA += 1 * digit(pdModes[i].channel, 10000);
		dA += 2 * digit(pdModes[i].channel, 1000);
		dA += 3 * digit(pdModes[i].channel, 100);
		dA += 3 * digit(pdModes[i].channel, 10);
		dA += 4 * digit(pdModes[i].channel, 1);

		double rate = interpolateEquidistant(lg, lgmin, lgmax, pdModes[i].rate);
		lossRate += rate * dA / A;
	}

	// comological scaling of interaction distance (in physical units)
	lossRate *= pow(1 + z, 3) * photonFieldScaling(photonField, z);
	return 1. / lossRate;
}

} // namespace crpropa

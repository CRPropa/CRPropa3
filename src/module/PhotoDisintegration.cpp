#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {


const double PhotoDisintegration::lgmin = 6.; // minimum log10(Lorentz-factor)
const double PhotoDisintegration::lgmax = 14.; // maximum log10(Lorentz-factor)
const size_t PhotoDisintegration::nlg = 201; // number of Lorentz-factor steps


PhotoDisintegration::PhotoDisintegration(PhotonField f, double limit) {
	setPhotonField(f);
	this->limit = limit;
}

void PhotoDisintegration::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("PhotoDisintegration: CMB");
		initRate(getDataPath("pd_CMB.txt"));
		initBranching(getDataPath("pd_branching_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("PhotoDisintegration: IRB (Kneiske 2004)");
		initRate(getDataPath("pd_IRB_Kneiske04.txt"));
		initBranching(getDataPath("pd_branching_IRB_Kneiske04.txt"));
		break;
	case IRB_Stecker05:
		setDescription("PhotoDisintegration: IRB (Stecker 2005)");
		initRate(getDataPath("pd_IRB_Stecker05.txt"));
		initBranching(getDataPath("pd_branching_IRB_Stecker05.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("PhotoDisintegration: IRB (Franceschini 2008)");
		initRate(getDataPath("pd_IRB_Franceschini08.txt"));
		initBranching(getDataPath("pd_branching_IRB_Franceschini08.txt"));
		break;
	case IRB_Finke10:
		setDescription("PhotoDisintegration: IRB (Finke 2010)");
		initRate(getDataPath("pd_IRB_Finke10.txt"));
		initBranching(getDataPath("pd_branching_IRB_Finke10.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("PhotoDisintegration: IRB (Dominguez 2011)");
		initRate(getDataPath("pd_IRB_Dominguez11.txt"));
		initBranching(getDataPath("pd_branching_IRB_Dominguez11.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("PhotoDisintegration: IRB (Gilmore 12)");
		initRate(getDataPath("pd_IRB_Gilmore12.txt"));
		initBranching(getDataPath("pd_branching_IRB_Gilmore12.txt"));
		break;
	default:
		throw std::runtime_error(
				"PhotoDisintegration: unknown photon background");
	}
}

void PhotoDisintegration::setLimit(double limit) {
	this->limit = limit;
}

void PhotoDisintegration::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error(
				"PhotoDisintegration: could not open file " + filename);

	// clear previously loaded interaction rates
	pdRate.clear();
	pdRate.resize(27 * 31);

	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;
		std::stringstream lineStream(line);

		int Z, N;
		lineStream >> Z;
		lineStream >> N;

		double r;
		for (size_t i = 0; i < nlg; i++) {
			lineStream >> r;
			pdRate[Z * 31 + N].push_back(r / Mpc);
		}
	}
	infile.close();
}

void PhotoDisintegration::initBranching(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error(
				"PhotoDisintegration: could not open file " + filename);

	// clear previously loaded interaction rates
	pdBranch.clear();
	pdBranch.resize(27 * 31);

	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;

		std::stringstream lineStream(line);

		int Z, N;
		lineStream >> Z;
		lineStream >> N;

		Branch branch;
		lineStream >> branch.channel;

		double r;
		for (size_t i = 0; i < nlg; i++) {
			lineStream >> r;
			branch.branchingRatio.push_back(r);
		}

		pdBranch[Z * 31 + N].push_back(branch);
	}

	infile.close();
}

void PhotoDisintegration::process(Candidate *candidate) const {
	// execute the loop at least once for limiting the next step
	double step = candidate->getCurrentStep();
	do {
		// check if nucleus
		int id = candidate->current.getId();
		if (not (isNucleus(id)))
			return;

		int A = massNumber(id);
		int Z = chargeNumber(id);
		int N = A - Z;
		size_t idx = Z * 31 + N;

		// check if disintegration data available
		if ((Z > 26) or (N > 30))
			return;
		if (pdRate[idx].size() == 0)
			return;

		// check if in tabulated energy range
		double z = candidate->getRedshift();
		double lg = log10(candidate->current.getLorentzFactor() * (1 + z));
		if ((lg <= lgmin) or (lg >= lgmax))
			return;

		double rate = interpolateEquidistant(lg, lgmin, lgmax, pdRate[idx]);

		// cosmological scaling, rate per comoving distance
		rate *= pow(1 + z, 2) * photonFieldScaling(photonField, z);

		Random &random = Random::instance();
		double randDistance = -log(random.rand()) / rate;

		// check if an interaction occurs in this step
		// if not, limit next step to a fraction of the mean free path
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}

		// index of closest tabulation point
		int l = round((lg - lgmin) / (lgmax - lgmin) * (nlg - 1));

		// select channel and interact
		const std::vector<Branch> &branches = pdBranch[idx];
		double cmp = random.rand();
		size_t i = 0;
		while ((i < branches.size()) and (cmp > 0)) {
			cmp -= branches[i].branchingRatio[l];
			i++;
		}
		performInteraction(candidate, branches[i-1].channel);

		// repeat with remaining step
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

double PhotoDisintegration::lossLength(int id, double gamma, double z) {
	// check if nucleus
	if (not (isNucleus(id)))
		return 0;

	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;
	size_t idx = Z * 31 + N;

	// check if disintegration data available
	if ((Z > 26) or (N > 30))
		return std::numeric_limits<double>::max();
	const std::vector<double> &rate = pdRate[idx];
	if (rate.size() == 0)
		return std::numeric_limits<double>::max();

	// check if in tabulated energy range
	double lg = log10(gamma) * (1 + z);
	if ((lg <= lgmin) or (lg >= lgmax))
		return std::numeric_limits<double>::max();

	// total interaction rate
	double lossRate = interpolateEquidistant(lg, lgmin, lgmax, rate);

	// comological scaling, rate per physical distance
	lossRate *= pow(1 + z, 3) * photonFieldScaling(photonField, z);

	// average number of nucleons lost for all disintegration channels
	double avg_dA = 0;
	const std::vector<Branch> &branches = pdBranch[idx];
	for (size_t i = 0; i < branches.size(); i++) {
		int channel = branches[i].channel;
		int dA = 0;
		dA += 1 * digit(channel, 100000);
		dA += 1 * digit(channel, 10000);
		dA += 2 * digit(channel, 1000);
		dA += 3 * digit(channel, 100);
		dA += 3 * digit(channel, 10);
		dA += 4 * digit(channel, 1);

		double br = interpolateEquidistant(lg, lgmin, lgmax,
				branches[i].branchingRatio);
		avg_dA += br * dA;
	}

	lossRate *= avg_dA / A;
	return 1 / lossRate;
}

} // namespace crpropa

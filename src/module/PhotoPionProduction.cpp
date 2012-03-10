#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/common.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace mpc {

PhotoPionProduction::PhotoPionProduction(PhotonField photonField) {
	init(photonField);
}

PhotoPionProduction::PhotoPionProduction() {
	init(CMBIR);
}

void PhotoPionProduction::init(PhotonField field) {
	name = "mpc::PhotoPionProduction";

	photonField = field;
	switch (photonField) {
	case CMB:
		init(getDataPath("PhotoPionProduction/cmb.txt"));
		break;
	case IR:
		init(getDataPath("PhotoPionProduction/ir.txt"));
		break;
	case CMBIR:
		init(getDataPath("PhotoPionProduction/cmbir.txt"));
		break;
	}
}

void PhotoPionProduction::init(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(name + ": could not open file " + filename);

	std::vector<double> x, yp, yn;
	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b, c;
			infile >> a >> b >> c;
			if (infile) {
				x.push_back(a * EeV); // energy in [EeV]
				yp.push_back(b / Mpc); // rate in [1/Mpc]
				yn.push_back(c / Mpc); // rate in [1/Mpc]
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	infile.close();

	acc = gsl_interp_accel_alloc();
	pRate = gsl_spline_alloc(gsl_interp_linear, x.size());
	nRate = gsl_spline_alloc(gsl_interp_linear, x.size());
	gsl_spline_init(pRate, &x[0], &yp[0], x.size());
	gsl_spline_init(nRate, &x[0], &yn[0], x.size());
}

PhotoPionProduction::~PhotoPionProduction() {
	gsl_spline_free(pRate);
	gsl_spline_free(nRate);
	gsl_interp_accel_free(acc);
}

std::string PhotoPionProduction::getDescription() const {
	switch (photonField) {
	case CMB: {
		return "Photo-pion production on CMB";
		break;
	}
	case IR: {
		return "Photo-pion production on IR";
		break;
	}
	case CMBIR: {
		return "Photo-pion production on CMBIR";
		break;
	}
	default: {
		return "Photo-pion production (unknown)";
		break;
	}
	}
}

void PhotoPionProduction::process(Candidate *candidate) const {
	double step = candidate->getCurrentStep();
	InteractionState interaction;

	while (true) {
		// check if interaction is set
		bool noState = !candidate->getInteractionState(name, interaction);
		if (noState) {
			// try to set a new interaction
			bool noInteraction = !setNextInteraction(candidate);
			if (noInteraction)
				return;
			// get the new interaction
			candidate->getInteractionState(name, interaction);
		}

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

bool PhotoPionProduction::setNextInteraction(Candidate *candidate) const {
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	int N = A - Z;
	double EpA = E / A * (1 + z); // CMB energies increase with (1+z)^3

	// check if out of energy range
	if ((EpA < 10 * EeV) or (EpA > 1e5 * EeV))
		return false;

	// find interaction with minimum random distance
	// check for interaction on protons
	InteractionState interaction;
	interaction.distance = std::numeric_limits<double>::max();
	Random &random = Random::instance();
	if (Z > 0) {
		double rate = gsl_spline_eval(pRate, EpA, acc) * Z;
		interaction.distance = -log(random.rand()) / rate;
		interaction.channel = 1;
	}
	// check for interaction on neutrons
	if (N > 0) {
		double rate = gsl_spline_eval(nRate, EpA, acc) * N;
		double d = -log(random.rand()) / rate;
		if (d < interaction.distance) {
			interaction.distance = -log(random.rand()) / rate;
			interaction.channel = 0;
		}
	}

	// CMB density increases with (1+z)^3 -> free distance decreases accordingly
	interaction.distance /= pow((1 + z), 3);

	if (isnan(interaction.distance)) {
		return false;
	}

	candidate->setInteractionState(name, interaction);
	return true;
}

void PhotoPionProduction::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(name, interaction);
	candidate->clearInteractionStates();

	// charge number loss of interaction nucleus
	int dZ = interaction.channel;
	// final proton number of emitted nucleon
	int Zfinal = dZ;
	// 1/3 probability of isospin change p <-> n
	if (Random::instance().rand() < 1. / 3.)
		Zfinal = abs(Zfinal - 1);

	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();

	if (A == 1) {
		// interaction on single nucleon
		candidate->current.setEnergy(E * 938. / 1232.);
		candidate->current.setId(getNucleusId(1, Zfinal));
	} else {
		// interaction on nucleus, update nucleus and emit nucleon
		candidate->current.setEnergy(E * (A - 1) / A);
		candidate->current.setId(getNucleusId(A - 1, Z - dZ));
		candidate->addSecondary(getNucleusId(1, Zfinal), E / A * 938. / 1232.);
	}
}

} // namespace mpc

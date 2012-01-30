#include "mpc/module/PhotoPionProduction.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>

namespace mpc {

PhotoPionProduction::PhotoPionProduction(PhotonField photonField) {
	init(photonField);
}

PhotoPionProduction::PhotoPionProduction() {
	init(CMBIR);
}

void PhotoPionProduction::init(PhotonField field) {
	photonField = field;
	switch (photonField) {
	case CMB:
		init("data/PhotoPionProduction/cmb.txt");
		break;
	case IR:
		init("data/PhotoPionProduction/ir.txt");
		break;
	case CMBIR:
		init("data/PhotoPionProduction/cmbir.txt");
		break;
	}
}

void PhotoPionProduction::init(std::string filename) {
	gsl_rng_env_setup();
	rand = gsl_rng_alloc(gsl_rng_default);

	// load interaction rate table
	std::vector<double> x, yp, yn;
	std::ifstream infile(filename.c_str());
	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b, c;
			infile >> a >> b >> c;
			if (infile) {
				x.push_back(a * EeV);
				yp.push_back(Mpc / b);
				yn.push_back(Mpc / c);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
	acc = gsl_interp_accel_alloc();
	protonFreePath = gsl_spline_alloc(gsl_interp_linear, x.size());
	neutronFreePath = gsl_spline_alloc(gsl_interp_linear, x.size());
	gsl_spline_init(protonFreePath, &x[0], &yp[0], x.size());
	gsl_spline_init(neutronFreePath, &x[0], &yn[0], x.size());
}

PhotoPionProduction::~PhotoPionProduction() {
	gsl_spline_free(protonFreePath);
	gsl_spline_free(neutronFreePath);
	gsl_interp_accel_free(acc);
	gsl_rng_free(rand);
}

std::string PhotoPionProduction::getDescription() const {
	switch (photonField) {
	case CMB: {
		return "Photo-pion production (CMB)";
		break;
	}
	case IR: {
		return "Photo-pion production (IR)";
		break;
	}
	case CMBIR: {
		return "Photo-pion production (CMBIR)";
		break;
	}
	default: {
		return "Photo-pion production (unknown)";
		break;
	}
	}
}

void PhotoPionProduction::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	int id = candidate->current.getId();

	// check if a interaction is already cached, if not set it
	if (id != cached_id) {
		cached_id = id;
		setNextInteraction(candidate);
	}

	// reduce distance to interaction
	cached_distance -= candidate->getCurrentStep();

	// check if free path is over
	if (cached_distance > 0) {
		candidate->limitNextStep(cached_distance);
		return;
	}
	performInteraction(candidate, secondaries);
	cached_id = 0;
	// should repeat this, goto?
}

void PhotoPionProduction::setNextInteraction(Candidate *candidate) {
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	int A = getMassNumberFromNucleusId(id);
	int Z = getChargeNumberFromNucleusId(id);
	int N = A - Z;

	double xi = E / A * (1 + z);
	double d, mfp;
	cached_distance = std::numeric_limits<double>::max();
	// check for interaction on proton
	if (Z > 0) {
		double mfp = gsl_spline_eval(protonFreePath, xi, acc) / Z;
		double d = -log(gsl_rng_uniform(rand)) * mfp;
		if (d < cached_distance) {
			cached_distance = d;
			cached_interaction = 1;
		}
	}
	// check for interaction on neutron
	if (N > 0) {
		double mfp = gsl_spline_eval(neutronFreePath, xi, acc) / N;
		double d = -log(gsl_rng_uniform(rand)) * mfp;
		if (d < cached_distance) {
			cached_distance = d;
			cached_interaction = 0;
		}
	}
}

void PhotoPionProduction::performInteraction(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	double E = candidate->current.getEnergy();
	int id = candidate->current.getId();
	int A = getMassNumberFromNucleusId(id);
	int Z = getChargeNumberFromNucleusId(id);

	int Zfinal = cached_interaction; // = 1, 0 for interaction on p, n
	if (gsl_rng_uniform(rand) > 2. / 3)
		Zfinal = abs(Zfinal - 1); // 2/3 probability of isospin change p <-> n

	// interaction on single nucleon
	if (A == 1) {
		candidate->current.setEnergy(E * 938. / 1232.);
		candidate->current.setId(getNucleusId(Zfinal, 1));
		return;
	}

	// else, interaction on nucleus
	candidate->current.setEnergy(E * (A - 1) / A);
	candidate->current.setId(getNucleusId(Z - cached_interaction, A - 1));
	candidate->addSecondary(getNucleusId(Zfinal, 1), E / A * 938. / 1232.);
}

} // namespace mpc

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
	std::vector<double> x, yp, yn;
	std::ifstream infile(filename.c_str());
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

void PhotoPionProduction::process(Candidate *candidate) {
	double step = candidate->getCurrentStep();

	while (true) {
		int id = candidate->current.getId();

		// set a new interaction if necessary
		if (id != cached_id) {
			// return if no data
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
		// counter over: interact
		cached_id = 0;
		step -= cached_distance;
		performInteraction(candidate);
	}
}

bool PhotoPionProduction::setNextInteraction(Candidate *candidate) {
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	int A = getMassNumberFromNucleusId(id);
	int Z = getChargeNumberFromNucleusId(id);
	int N = A - Z;
	double EpA = E / A * (1 + z); // CMB energies increase with (1+z)^3

	// out of energy range
	if ((EpA < 10 * EeV) or (EpA > 1e5 * EeV))
		return false;

	cached_distance = std::numeric_limits<double>::max();
	// check for interaction on proton
	if (Z > 0) {
		double rate = gsl_spline_eval(pRate, EpA, acc) * Z;
		cached_distance = -log(mtrand.rand()) / rate;
		cached_interaction = 1;
	}
	// check for interaction on neutron
	if (N > 0) {
		double rate = gsl_spline_eval(nRate, EpA, acc) * N;
		double d = -log(mtrand.rand()) / rate;
		if (d < cached_distance) {
			cached_distance = d;
			cached_interaction = 0;
		}
	}

	cached_distance /= pow((1 + z), 3); // CMB density increases with (1+z)^3
	return true;
}

void PhotoPionProduction::performInteraction(Candidate *candidate) {
	double E = candidate->current.getEnergy();
	int id = candidate->current.getId();
	int A = getMassNumberFromNucleusId(id);
	int Z = getChargeNumberFromNucleusId(id);

	// final proton number of interaction nucleon
	int Zfinal = cached_interaction;
	if (mtrand.rand() < 1. / 3.)
		Zfinal = abs(Zfinal - 1); // 1/3 probability of isospin change p <-> n

	// interaction on single nucleon
	if (A == 1) {
		candidate->current.setEnergy(E * 938. / 1232.);
		candidate->current.setId(getNucleusId(1, Zfinal));
		return;
	}

	// else, interaction on nucleus, nucleon is emitted
	candidate->current.setEnergy(E * (A - 1) / A);
	candidate->current.setId(getNucleusId(A - 1, Z - cached_interaction));
	candidate->addSecondary(getNucleusId(Zfinal, 1), E / A * 938. / 1232.);
}

} // namespace mpc

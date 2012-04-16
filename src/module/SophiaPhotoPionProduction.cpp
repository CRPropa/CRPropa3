#include "mpc/module/SophiaPhotoPionProduction.h"

#include "sophia.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace mpc {

SophiaPhotoPionProduction::SophiaPhotoPionProduction(PhotonField photonField) {
	init(photonField);
}

SophiaPhotoPionProduction::SophiaPhotoPionProduction() {
	init(CMB);
}

void SophiaPhotoPionProduction::init(PhotonField field) {
	name = "mpc::PhotoPionProduction";

	photonField = field;
	switch (photonField) {
	case CMB:
		init(getDataPath("PhotoPionProduction/cmb.txt"));
		break;
	case IR:
		init(getDataPath("PhotoPionProduction/ir.txt"));
		break;
	}
}

void SophiaPhotoPionProduction::init(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(name + ": could not open file " + filename);

	std::vector<double> x, yp, yn;
	x.reserve(100);
	yp.reserve(100);
	yn.reserve(100);
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
	Emin = x.front();
	Emax = x.back();
}

SophiaPhotoPionProduction::~SophiaPhotoPionProduction() {
	gsl_spline_free(pRate);
	gsl_spline_free(nRate);
	gsl_interp_accel_free(acc);
}

std::string SophiaPhotoPionProduction::getDescription() const {
	switch (photonField) {
	case CMB: {
		return "Photo-pion production on CMB using SOPHIA";
		break;
	}
	case IR: {
		return "Photo-pion production on IR using SOPHIA";
		break;
	}
	default: {
		return "Photo-pion production using SOPHIA (unknown)";
		break;
	}
	}
}

void SophiaPhotoPionProduction::process(Candidate *candidate) const {
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

bool SophiaPhotoPionProduction::setNextInteraction(Candidate *candidate) const {
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	int N = A - Z;
	double EpA = E / A * (1 + z); // CMB energies increase with (1+z)^3

	// check if out of energy range
	if ((EpA < Emin) or (EpA > Emax))
		return false;

	// find interaction with minimum random distance
	// check for interaction on protons
	InteractionState interaction;
	interaction.distance = std::numeric_limits<double>::max();
	Random &random = Random::instance();
	if (Z > 0) {
		double rate = gsl_spline_eval(pRate, EpA, acc) * Z;
		if (rate != 0) {
			interaction.distance = -log(random.rand()) / rate;
			interaction.channel = 1;
		}
	}
	// check for interaction on neutrons
	if (N > 0) {
		double rate = gsl_spline_eval(nRate, EpA, acc) * N;
		if (rate != 0) {
			double d = -log(random.rand()) / rate;
			if (d < interaction.distance) {
				interaction.distance = -log(random.rand()) / rate;
				interaction.channel = 0;
			}
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

void SophiaPhotoPionProduction::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(name, interaction);
	candidate->clearInteractionStates();

	int dZ = interaction.channel; // charge number loss of interacting nucleus

	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	double EpA = E / A;

	// arguments for sophia
	int nature = 1 - dZ; // interacting particle: 0 for proton, 1 for neutron
	double Ein = EpA / GeV; // energy of in-going nucleon in GeV
	double momentaList[5][2000]; // momentum list, what are the five components?
	int particleList[2000]; // particle id list
	int nParticles; // number of outgoing particles
	double redshift = candidate->getRedshift();
	int background = photonField + 1; // Photon background: 1 for CMB, 2 for Primack IR
	double maxRedshift = 100; // IR photon density is zero above this redshift
	int IRNE = 2; // ???
	double dummy[2]; // ???

	sophiaevent_(nature, Ein, momentaList, particleList, nParticles, redshift,
			background, maxRedshift, IRNE, dummy, dummy);

	for (int i = 0; i < nParticles; i++) { // loop over out-going particles
		double Eout = momentaList[3][i] * GeV; // only the energy is used; could be changed for more realism
		int pType = particleList[i];
//		13  proton
//		14  neutron
//		-13 antiproton
//		-14 antineutron
//		1   photon
//		2   e+
//		3   e-
//		15  nu_e
//		16  antinu_e
//		17  nu_muon
//		18  antinu_muon

		if (pType == 13 or pType == 14) {
			int Zout = 14 - pType; // 1 for proton, 0 for neutron
			if (A == 1) { // in-going particle was a nucleon: update its properties
				candidate->current.setEnergy(Eout);
				candidate->current.setId(getNucleusId(1, Zout));
			} else {
				// // in-going particle was a nucleus: update nucleus and emit nucleon
				candidate->current.setEnergy(E - Eout);
				candidate->current.setId(getNucleusId(A - 1, Z - dZ));
				candidate->addSecondary(getNucleusId(1, Zout), Eout);
			}
		} else if (pType == -13 or pType == -14) {
			continue; // anti-proton/neutron annihilate quickly -> new electromagnetic cascade?
		} else if (pType >= 1 and pType <= 3) {
			continue; // new electromagnetic cascade
		} else if (pType >= 15) {
			continue; // new neutrino
		}
	}
}

} // namespace mpc

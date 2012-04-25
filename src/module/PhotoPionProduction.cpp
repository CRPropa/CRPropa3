#include "mpc/module/PhotoPionProduction.h"

#include <kiss/convert.h>
#include "sophia.h"

#include <limits>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace mpc {

PhotoPionProduction::PhotoPionProduction(int photonField) {
	init(photonField);
}

void PhotoPionProduction::init(int photonField) {
	this->photonField = photonField;
	switch (photonField) {
	case CMB:
		setDescription("PhotoPionProduction:CMB");
		init(getDataPath("PhotoPionProduction/PPtable_CMB.txt"));
		break;
	case IRB:
		setDescription("PhotoPionProduction:IRB");
		init(getDataPath("PhotoPionProduction/PPtable_IRB.txt"));
		break;
	case CMB_IRB:
		setDescription("PhotoPionProduction:CMB_IRB");
		init(getDataPath("PhotoPionProduction/PPtable_CMB_IRB.txt"));
		break;
	default:
		throw std::runtime_error(
				"mpc::PhotoPionProduction: unknown photon background");
	}
}

void PhotoPionProduction::init(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"mpc::PhotoPionProduction: could not open file " + filename);

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

PhotoPionProduction::~PhotoPionProduction() {
	gsl_spline_free(pRate);
	gsl_spline_free(nRate);
	gsl_interp_accel_free(acc);
}

bool PhotoPionProduction::setNextInteraction(Candidate *candidate,
		InteractionState &interaction) const {
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

	candidate->setInteractionState(getDescription(), interaction);
	return true;
}

void PhotoPionProduction::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(getDescription(), interaction);
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

SophiaPhotoPionProduction::SophiaPhotoPionProduction(int photonField) {
	if ((photonField == CMB) or (photonField == IRB))
		init(photonField);
	else
		std::runtime_error("mpc::SophiaPhotoPionProduction: only CMB or IRB background possible");
}

void SophiaPhotoPionProduction::performInteraction(Candidate *candidate) const {
	InteractionState interaction;
	candidate->getInteractionState(getDescription(), interaction);
	candidate->clearInteractionStates();

	int channel = interaction.channel; // 1 for interaction proton, 0 for neutron

	double E = candidate->current.getEnergy();
	int A = candidate->current.getMassNumber();
	int Z = candidate->current.getChargeNumber();
	double EpA = E / A;

	// arguments for sophia
	int nature = 1 - channel; // interacting particle: 0 for proton, 1 for neutron
	double Ein = EpA / GeV; // energy of in-going nucleon in GeV
	double momentaList[5][2000]; // momentum list, what are the five components?
	int particleList[2000]; // particle id list
	int nParticles; // number of outgoing particles
	double redshift = candidate->getRedshift();
	int background; // Photon background: 1 for CMB, 2 for Kneiske IRB
	switch (photonField) {
	case CMB:
		background = 1;
		break;
	case IRB:
		background = 2;
		break;
	}
	double maxRedshift = 100; // IR photon density is zero above this redshift
	int dummy1;
	double dummy2[2];

	sophiaevent_(nature, Ein, momentaList, particleList, nParticles, redshift,
			background, maxRedshift, dummy1, dummy2, dummy2);

	for (int i = 0; i < nParticles; i++) { // loop over out-going particles
		double Eout = momentaList[3][i] * GeV; // only the energy is used; could be changed for more detail
		int pType = particleList[i];

		switch (pType) {
		case 13: // proton
		case 14: // neutron
			if (A == 1) { // in-going particle was a nucleon: update its properties
				candidate->current.setEnergy(Eout);
				candidate->current.setId(getNucleusId(1, 14 - pType));
			} else { // in-going particle was a nucleus: update nucleus and emit nucleon
				candidate->current.setEnergy(E - Eout);
				candidate->current.setId(getNucleusId(A - 1, Z - channel));
				candidate->addSecondary(getNucleusId(1, 14 - pType), Eout);
			}
			break;
		case -13: // anti-proton
			candidate->addSecondary(-getNucleusId(1, 1), Eout);
			break;
		case -14: // anti-neutron
			candidate->addSecondary(-getNucleusId(1, 0), Eout);
			break;
		case 1: // photon
			candidate->addSecondary(22, Eout);
			break;
		case 2: // positron
			candidate->addSecondary(-11, Eout);
			break;
		case 3: // electron
			candidate->addSecondary(11, Eout);
			break;
		case 15: // nu_e
			candidate->addSecondary(12, Eout);
			break;
		case 16: // antinu_e
			candidate->addSecondary(-12, Eout);
			break;
		case 17: // nu_muon
			candidate->addSecondary(14, Eout);
			break;
		case 18: // antinu_muon
			candidate->addSecondary(-14, Eout);
			break;
		default:
			throw std::runtime_error(
					"mpc::PhotoPionProduction: unexpected particle "
							+ kiss::str(pType));
		}
	}
}

} // namespace mpc

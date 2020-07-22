#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

ElectronPairProduction::ElectronPairProduction(PhotonField photonField,
		bool haveElectrons, double limit) {
	setPhotonField(photonField);
	this->haveElectrons = haveElectrons;
	this->limit = limit;
}

void ElectronPairProduction::setPhotonField(PhotonField photonField) {
	this->photonField = photonField;
	std::string fname = photonField.getFieldName();
	setDescription("ElectronPairProduction: " + fname);
	initRate(getDataPath("ElectronPairProduction/lossrate_" + fname + ".txt"));
	initSpectrum(getDataPath("ElectronPairProduction/spectrum_" + fname.substr(0,3) + ".txt"));
}

void ElectronPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void ElectronPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void ElectronPairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("ElectronPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabLorentzFactor.clear();
	tabLossRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabLorentzFactor.push_back(pow(10, a));
				tabLossRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void ElectronPairProduction::initSpectrum(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("ElectronPairProduction: could not open file " + filename);

	double dNdE;
	tabSpectrum.resize(70);
	for (size_t i = 0; i < 70; i++) {
		tabSpectrum[i].resize(170);
		for (size_t j = 0; j < 170; j++) {
			infile >> dNdE;
			tabSpectrum[i][j] = dNdE * pow(10, (7 + 0.1 * j)); // read electron distribution pdf(Ee) ~ dN/dEe * Ee
		}
		for (size_t j = 1; j < 170; j++) {
			tabSpectrum[i][j] += tabSpectrum[i][j - 1]; // cdf(Ee), unnormalized
		}
	}
	infile.close();
}

double ElectronPairProduction::lossLength(int id, double lf, double z) const {
	double Z = chargeNumber(id);
	if (Z == 0)
		return std::numeric_limits<double>::max(); // no pair production on uncharged particles

	lf *= (1 + z);
	if (lf < tabLorentzFactor.front())
		return std::numeric_limits<double>::max(); // below energy threshold

	double rate;
	if (lf < tabLorentzFactor.back())
		rate = interpolate(lf, tabLorentzFactor, tabLossRate); // interpolation
	else
		rate = tabLossRate.back() * pow(lf / tabLorentzFactor.back(), -0.6); // extrapolation

	double A = nuclearMass(id) / mass_proton; // more accurate than massNumber(Id)
	rate *= Z * Z / A * pow(1 + z, 3) * photonField.getRedshiftScaling(z);
	return 1. / rate;
}

void ElectronPairProduction::process(Candidate *c) const {
	int id = c->current.getId();
	if (not (isNucleus(id)))
		return; // only nuclei

	double lf = c->current.getLorentzFactor();
	double z = c->getRedshift();
	double losslen = lossLength(id, lf, z);  // energy loss length
	if (losslen >= std::numeric_limits<double>::max())
		return;

	double step = c->getCurrentStep() / (1 + z); // step size in local frame
	double loss = step / losslen;  // relative energy loss

	if (haveElectrons) {
		double dE = c->current.getEnergy() * loss;  // energy loss
		int i = round((log10(lf) - 6.05) * 10);  // find closest cdf(Ee|log10(gamma))
		i = std::min(std::max(i, 0), 69);
		Random &random = Random::instance();

		// draw pairs as long as their energy is smaller than the pair production energy loss
		while (dE > 0) {
			size_t j = random.randBin(tabSpectrum[i]);
			double Ee = pow(10, 6.95 + (j + random.rand()) * 0.1) * eV;
			double Epair = 2 * Ee; // NOTE: electron and positron in general don't have same lab frame energy, but averaged over many draws the result is consistent
			// if the remaining energy is not sufficient check for random accepting
			if (Epair > dE)
				if (random.rand() > (dE / Epair))
					break; // not accepted

			// create pair and repeat with remaining energy
			dE -= Epair;
			Vector3d pos = random.randomInterpolatedPosition(c->previous.getPosition(), c->current.getPosition());
			c->addSecondary( 11, Ee, pos);
			c->addSecondary(-11, Ee, pos);
		}
	}

	c->current.setLorentzFactor(lf * (1 - loss));
	c->limitNextStep(limit * losslen);
}

} // namespace crpropa

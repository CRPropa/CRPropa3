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
	switch (photonField) {
	case CMB:
		setDescription("ElectronPairProduction: CMB");
		initRate(getDataPath("epp_CMB.txt"));
		initSpectrum(getDataPath("epp_spectrum_CMB.txt"));
		break;
	case IRB:  // default: Kneiske '04 IRB model
	case IRB_Kneiske04:
		setDescription("ElectronPairProduction: IRB (Kneiske 2004)");
		initRate(getDataPath("epp_IRB_Kneiske04.txt"));
		initSpectrum(getDataPath("epp_spectrum_IRB.txt"));
		break;
	case IRB_Stecker05:
		setDescription("ElectronPairProduction: IRB (Stecker 2005)");
		initRate(getDataPath("epp_IRB_Stecker05.txt"));
		initSpectrum(getDataPath("epp_spectrum_IRB.txt"));
		break;
	case IRB_Franceschini08:
		setDescription("ElectronPairProduction: IRB (Franceschini 2008)");
		initRate(getDataPath("epp_IRB_Franceschini08.txt"));
		initSpectrum(getDataPath("epp_spectrum_IRB.txt"));
		break;
	case IRB_Finke10:
		setDescription("ElectronPairProduction: IRB (Finke 2010)");
		initRate(getDataPath("epp_IRB_Finke10.txt"));
		initSpectrum(getDataPath("epp_spectrum_IRB.txt"));
		break;
	case IRB_Dominguez11:
		setDescription("ElectronPairProduction: IRB (Dominguez 2011)");
		initRate(getDataPath("epp_IRB_Dominguez11.txt"));
		initSpectrum(getDataPath("epp_spectrum_IRB.txt"));
		break;
	case IRB_Gilmore12:
		setDescription("ElectronPairProduction: IRB (Gilmore 2012)");
		initRate(getDataPath("epp_IRB_Gilmore12.txt"));
		initSpectrum(getDataPath("epp_spectrum_IRB.txt"));
		break;
	default:
		throw std::runtime_error(
				"ElectronPairProduction: unknown photon background");
	}
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
		throw std::runtime_error(
				"ElectronPairProduction: could not open file " + filename);

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
	// proton energies: 70 energies from 10^15 - 10^22 eV
	for (size_t i = 0; i < 70; i++)
		tabE.push_back(pow(10, 15 + i * 7. / 69.) * eV);
	// electron energies: 171 energies from 10^6.95 - 10^23.95 eV
	for (size_t j = 0; j < 171; j++)
		tabEe.push_back(pow(10, 6.95 + 0.1 * j) * eV);
	for (size_t j = 0; j < 170; j++)
		tabEeWidth.push_back(tabEe[j+1] - tabEe[j]);

	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error(
				"ElectronPairProduction: could not open file " + filename);

	double dNdE;
	tabSpectrum.resize(70);
	for (size_t i = 0; i < 70; i++) {
		tabSpectrum[i].resize(170);
		for (size_t j = 0; j < 170; j++) {
			infile >> dNdE;
			tabSpectrum[i][j] = dNdE * tabEe[j]; // read in spectrum: f(E) ~ dN/dE * E
		}
		for (size_t j = 1; j < 170; j++) {
			tabSpectrum[i][j] += tabSpectrum[i][j - 1]; // cumulate F(E), this does not need to be normalized
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
	rate *= Z * Z / A * pow(1 + z, 3) * photonFieldScaling(photonField, z);
	return 1. / rate;
}

void ElectronPairProduction::addElectrons(Candidate *c, double loss) const {
	double E = c->current.getEnergy();
	double dE = E * loss; // energy loss
	double z = c->getRedshift();
	double Eeff = E / massNumber(c->current.getId()) * (1 + z);

	// interpolate spectrum in the Eff
	size_t i = std::upper_bound(tabE.begin(), tabE.end(), Eeff) - tabE.begin() - 1;
	double a = (Eeff - tabE[i]) / (tabE[i + 1] - tabE[i]);

	std::vector<double> spectrum(170);
	for (size_t j = 0; j < 170; j++)
		spectrum[j] = tabSpectrum[i][j]
				+ a * (tabSpectrum[i + 1][j] - tabSpectrum[i][j]);

	// draw pairs as long as their energy is smaller than the pair production energy loss
	Random &random = Random::instance();
	while (dE > 0) {
		size_t i = random.randBin(spectrum); // draw random bin
		double Ee = tabEe[i] + random.rand() * tabEeWidth[i]; // draw random uniform energy in bin
		Ee /= (1 + z); // dN/dE(Ep,Ee,z) = (1+z)^4 * dN/dE(Ep*(1+z),Ee*(1+z),0)

		double Epair = 2 * Ee; // electron and positron should generally not have the same energy in the lab frame, however averaged over many draws the result of this method will be consistent with the correct implementation

		// if the remaining energy is not sufficient check for random accepting
		if (Epair > dE)
			if (random.rand() > (dE / Epair))
				break; // not accepted

		// create pair and repeat with remaining energy
		dE -= Epair;
		c->addSecondary(11, Ee);
		c->addSecondary(-11, Ee);
	}
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

	if (haveElectrons)
		addElectrons(c, loss);

	c->current.setLorentzFactor(lf * (1 - loss));
	c->limitNextStep(limit * losslen);
}

} // namespace crpropa

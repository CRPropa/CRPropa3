#include "crpropa/module/PhotonDINT1D.h"
#include "crpropa/Units.h"
#include "crpropa/Cosmology.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "dint/prop_second.h"

using namespace std;

namespace crpropa {

class PhotonDINT1DImpl {
public:
	PhotonDINT1DImpl() {
		// Initialize the energy grids for dint
		New_dCVector(&energyGrid, NUM_MAIN_BINS);
		New_dCVector(&energyWidth, NUM_MAIN_BINS);
		SetEnergyBins(MIN_ENERGY_EXP, &energyGrid, &energyWidth);
	}

	virtual ~PhotonDINT1DImpl() {
		Delete_dCVector(&energyGrid);
		Delete_dCVector(&energyWidth);
	}

	dCVector energyGrid, energyWidth;
	virtual void saveSpectrum(Spectrum *spectrum) = 0;
};

class PhotonDINT1DAsciiImpl: public PhotonDINT1DImpl {
	mutable std::ofstream fout;
public:
	PhotonDINT1DAsciiImpl(const std::string &filename) :
			PhotonDINT1DImpl(), fout(filename.c_str()) {
		for (int j = 0; j < energyGrid.dimension; j++) {
			fout << (energyGrid.vector[j] * ELECTRON_MASS) << " ";
		}
		fout << std::endl;
	}

	void saveSpectrum(Spectrum *spectrum) {
		for (int j = 0; j < spectrum->numberOfMainBins; j++) {
			fout << spectrum->spectrum[PHOTON][j] << " "; // spectrum: mean number of particles per energy bin
		}
		fout << endl;
	}
};

#ifdef CRPROPA_HAVE_ROOT
class PhotonDINT1DROOTImpl: public PhotonDINT1DImpl {
public:
	PhotonDINT1DROOTImpl(const std::string &filename) :
	PhotonDINT1DImpl() {
		//TODO: new histogram
	}

	~PhotonDINT1DROOTImpl() {
		//TODO: delete histogram
	}

	void saveSpectrum(Spectrum *spectrum) {
		//TODO: add to histogram
	}
};
#endif

PhotonDINT1D::PhotonDINT1D(const string &filename) :
		filename(filename), IRFlag(2), RadioFlag(2), Zmax(5), Cutcascade_Magfield(
				0), impl(0) {
	dataPath = getDataPath("dint");

	std::string::size_type i = filename.find_last_of('.');
	std::string ext = filename.substr(i);
	if (ext == "root") {
#ifdef CRPROPA_HAVE_ROOT
		impl = new PhotonDINT1DROOTImpl(filename);
#else
		throw std::runtime_error("ROOT not compiled!");
#endif
	} else {
		impl = new PhotonDINT1DAsciiImpl(filename);
	}

}

PhotonDINT1D::~PhotonDINT1D() {
	delete impl;
}

void PhotonDINT1D::setIRFlag(int ir) {
	IRFlag = ir;
}

void PhotonDINT1D::setRadioFlag(int radio) {
	RadioFlag = radio;
}

void PhotonDINT1D::setZmax(double zmax) {
	Zmax = zmax;
}

void PhotonDINT1D::process(Candidate *candidate) const {
	if (candidate->current.getId() != 22)
		return;

	// Initialize the spectrum
	Spectrum inputSpectrum;
	NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);

	double criticalEnergy = candidate->current.getEnergy()
			/ (eV * ELECTRON_MASS); // units of dint
	int maxBin = (int) ((log10(criticalEnergy * ELECTRON_MASS) - MAX_ENERGY_EXP)
			* BINS_PER_DECADE + NUM_MAIN_BINS);
	inputSpectrum.spectrum[PHOTON][maxBin] = 1.;

	// Initialize the bField
	dCVector bField;
	New_dCVector(&bField, 1);

	// Initialize output spectrum
	Spectrum outputSpectrum;
	NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);

	double h = H0() * Mpc / 1000;
	double ol =  omegaL();
	double om = omegaM();
	double showerPropDistance = candidate->current.getPosition().getR() / Mpc;
	double z = candidate->getRedshift();
	if (z == 0) {
		//TODO: use z value for distance calculation
	}

	prop_second(showerPropDistance, &bField, &impl->energyGrid,
		    &impl->energyWidth, &inputSpectrum, &outputSpectrum, dataPath,
		    IRFlag, Zmax, RadioFlag, h, om, ol,
		    Cutcascade_Magfield);

#pragma omp critical
	{
		impl->saveSpectrum(&outputSpectrum);
	}

	DeleteSpectrum(&outputSpectrum);
	DeleteSpectrum(&inputSpectrum);
	Delete_dCVector(&bField);

	candidate->setActive(false);
}

string PhotonDINT1D::getDescription() const {
	std::stringstream s;
	s << "PhotonDINT1D: Output file = " << filename;
	return s.str();
}

} // namespace crpropa

#ifndef CRPROPA_EMINVERSECOMPTONSCATTERING_H
#define CRPROPA_EMINVERSECOMPTONSCATTERING_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include <fstream>

namespace crpropa {

/**
 @class EMInverseComptonScattering
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering of electrons with background photons for several photon fields.
 The upscattered photons are optionally created as secondary particles (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
*/
class EMInverseComptonScattering: public Module {
private:
	PhotonField photonField;
	ScalarGrid4d spaceTimeGrid;
	ScalarGrid spaceGrid;
	bool havePhotons;
	std::string tag;
	double limit;

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate

public:
	EMInverseComptonScattering(
		PhotonField photonField,	   //!< target photon background
		bool havePhotons = false,      //!< switch to create secondary photon
		std::string tag = "EMIC",
		double limit = 0.1             //!< step size limit as fraction of mean free path
	);

	EMInverseComptonScattering(
		PhotonField photonField,	   //!< target photon background
		ScalarGrid4d spaceTimeGrid,
		bool havePhotons = false,      //!< switch to create secondary photon
		std::string tag = "EMIC",
		double limit = 0.1             //!< step size limit as fraction of mean free path
	);

	EMInverseComptonScattering(
		PhotonField photonField,	   //!< target photon background
		ScalarGrid spaceGrid,
		bool havePhotons = false,      //!< switch to create secondary photon
		std::string tag = "EMIC",
		double limit = 0.1             //!< step size limit as fraction of mean free path
	);

	void setPhotonField(PhotonField photonField);
	void setHavePhotons(bool havePhotons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_EMINVERSECOMPTONSCATTERING_H

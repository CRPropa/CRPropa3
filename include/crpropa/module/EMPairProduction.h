#ifndef CRPROPA_EMPAIRPRODUCTION_H
#define CRPROPA_EMPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include <fstream>

namespace crpropa {

/**
 @class EMPairProduction
 @brief Electron-pair production of photons with background photons.

 This module simulates electron-pair production of cosmic ray photons with background photons:
 gamma + gamma_b --> e+ + e- (Breit-Wheeler process).
 The resulting electron positron pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 */
class EMPairProduction: public Module {
private:
	PhotonField photonField;
	ScalarGrid4d geometryGrid;
	bool haveElectrons;
	double limit;

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate

public:
	EMPairProduction(
		PhotonField photonField = CMB, //!< target photon background
		ScalarGrid4d geometryGrid = ScalarGrid4d(Vector3d(0.),0., 1,1,1,1, Vector3d(1.),1.),
		bool haveElectrons = false,    //!< switch to create secondary electron pair
		double limit = 0.1             //!< step size limit as fraction of mean free path
		);

	void setPhotonField(PhotonField photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void performInteraction(Candidate *candidate) const;
	void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_EMPAIRPRODUCTION_H

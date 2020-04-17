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
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class EMPairProduction: public Module {
private:
	PhotonField photonField;
	bool haveElectrons;
	double limit;
	double thinning;

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
		bool haveElectrons = false,    //!< switch to create secondary electron pair
		double thinning = 0,           //!< weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
		double limit = 0.1             //!< step size limit as fraction of mean free path
		);

	void setPhotonField(PhotonField photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);
	void setThinning(double thinning);

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void performInteraction(Candidate *candidate) const;
	void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_EMPAIRPRODUCTION_H

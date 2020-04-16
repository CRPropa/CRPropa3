#ifndef CRPROPA_EMDOUBLEPAIRPRODUCTION_H
#define CRPROPA_EMDOUBLEPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include <fstream>

namespace crpropa {

/**
 @class EMDoublePairProduction
 @brief Electron double pair production of photons with background photons.

 This module simulates electron double pair production of photons with background photons for several photon fields.
 The secondary electrons from this interaction are optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 */
class EMDoublePairProduction: public Module {
private:
	PhotonField photonField;
	bool haveElectrons;
	double limit;
	double thinning;

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]

public:
	EMDoublePairProduction(
		PhotonField photonField = CMB, //!< target photon background
		bool haveElectrons = false,    //!< switch to create the secondary electron pair
		double thinning = 0,           //!< weighted sampling of secondaries
		double limit = 0.1             //!< step size limit as fraction of mean free path
		);

	void setPhotonField(PhotonField photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);
	void setThinning(double thinning);

	void initRate(std::string filename);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};

} // namespace crpropa

#endif // CRPROPA_EMDOUBLEPAIRPRODUCTION_H

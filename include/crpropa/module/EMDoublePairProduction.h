#ifndef CRPROPA_EMDOUBLEPAIRPRODUCTION_H
#define CRPROPA_EMDOUBLEPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include <fstream>

namespace crpropa {

/**
 @class EMDoublePairProduction
 @brief Electron double pair production of photons with background photons.

 This module simulates electron double pair production in photon background photon interactions.\n
 Several photon fields can be selected.\n
 By default, the module limits the step size to 10% of the energy loss length of the particle.
 */
class EMDoublePairProduction: public Module {
private:
	PhotonField photonField;

	std::vector<double> tabInteractionRate; /*< tabulated interaction rate in [1/m] */
	std::vector<double> tabPhotonEnergy; /*< tabulated photon energy in [J] */
	double limit; ///< fraction of energy loss length to limit the next step
	bool haveElectrons;

public:
	EMDoublePairProduction(PhotonField photonField = CMB, bool haveElectrons =
			false, double limit = 0.1);

	void setPhotonField(PhotonField photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};

} // namespace crpropa

#endif // CRPROPA_EMDOUBLEPAIRPRODUCTION_H

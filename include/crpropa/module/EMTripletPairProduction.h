#ifndef CRPROPA_EMTRIPLETPAIRPRODUCTION_H
#define CRPROPA_EMTRIPLETPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include <fstream>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class EMTripletPairProduction
 @brief Electron triplet pair production of electrons with background photons.

 This module simulates electron triplet pair production of electrons with background photons for several photon fields.
 The secondary electrons from this interaction are optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
*/
class EMTripletPairProduction: public Module {
private:
	PhotonField photonField;
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
	EMTripletPairProduction(
		PhotonField photonField = CMB, //!< target photon background
		bool haveElectrons = false,    //!< switch to create secondary electron pair
		double limit = 0.1             //!< step size limit as fraction of mean free path
		);

	void setPhotonField(PhotonField photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMTRIPLETPAIRPRODUCTION_H

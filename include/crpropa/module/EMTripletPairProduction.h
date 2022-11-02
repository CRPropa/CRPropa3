#ifndef CRPROPA_EMTRIPLETPAIRPRODUCTION_H
#define CRPROPA_EMTRIPLETPAIRPRODUCTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

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
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
*/
class EMTripletPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool haveElectrons;
	double limit;
	double thinning;
	std::string interactionTag = "EMTP";

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate

public:
	/** Constructor
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 */
	EMTripletPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1);

	void setPhotonField(ref_ptr<PhotonField> photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);
	void setThinning(double thinning);

	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;
	
	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMTRIPLETPAIRPRODUCTION_H

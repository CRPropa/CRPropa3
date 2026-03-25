#ifndef CRPROPA_EMDOUBLEPAIRPRODUCTION_H
#define CRPROPA_EMDOUBLEPAIRPRODUCTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/InteractionRates.h"
#include "crpropa/Geometry.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class EMDoublePairProduction
 @brief Electron double pair production of photons with background photons.

 This module simulates electron double pair production of photons with background photons for several photon fields.
 The secondary electrons from this interaction are optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 The surface is defined to include the nodes of the grid contained within.
 */
class EMDoublePairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool haveElectrons;
	double limit;
	double thinning;
	ref_ptr<Surface> surface;
	std::string interactionTag = "EMDP";
	ref_ptr<InteractionRates> interactionRates;
public:
	/** Constructor
	 The object used to load, store and access to the interaction rates of the process is the interactionRates pointer.
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 @param surface     suface to enclose the grid nodes to be loaded
	 */
	EMDoublePairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1, ref_ptr<Surface> surface = nullptr);
	
	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);
	
	// decide if secondary electrons are added to the simulation
	void setHaveElectrons(bool haveElectrons);
	
	/** Limit the propagation step to a fraction of the mean free path
	 * @param limit fraction of the mean free path
	 */
	void setLimit(double limit);
	
	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);
	
	/** Apply a surface that confine the position dependent photon field region.
	 * @param surface closed surface to confine the grid to be  uploaded
	 */
	void setSurface(ref_ptr<Surface> surface);
	ref_ptr<Surface> getSurface() const;

	/** set a custom interaction rate
	 * With this function you can change the type of interaction rate,
	 * if you would for example like to change from a homogeneous to a position
	 * dependent interaction rate.
	 * @param intRates ref_ptr to a InteractionRates class
	 */
	void setInteractionRates(ref_ptr<InteractionRates> intRates);
	ref_ptr<InteractionRates> getInteractionRates() const;
	
	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	/** Loads the interaction rate
	 * This function loads the interaction rate from a given file/folder.
	 * @param path The name of the file/folder containing the interaction rates
	 */
	void initRate(std::string path);
	
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;
	
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMDOUBLEPAIRPRODUCTION_H

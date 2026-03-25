#ifndef CRPROPA_EMINVERSECOMPTONSCATTERING_H
#define CRPROPA_EMINVERSECOMPTONSCATTERING_H

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
 @class EMInverseComptonScattering
 @brief Inverse Compton scattering of electrons with background photons.

 This module simulates inverse Compton scattering of electrons with background photons for several photon fields.
 The upscattered photons are optionally created as secondary particles (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 The surface is defined to include the nodes of the grid contained within.
*/
class EMInverseComptonScattering: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool havePhotons;
	double limit;
	double thinning;
	ref_ptr<Surface> surface;
	std::string interactionTag = "EMIC";
	ref_ptr<InteractionRates> interactionRates;	
public:
	/** Constructor
	 The object used to load, store and access to the interaction rates of the process is the interactionRates pointer.
	 @param photonField		target photon field
	 @param havePhotons		if true, add secondary photons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 @param surface			suface to enclose the grid nodes to be loaded
	 */
	EMInverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons = false, double thinning = 0, double limit = 0.1, ref_ptr<Surface> surface = nullptr);
	
	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);
	
	// decide if secondary photons are added to the simulation
	void setHavePhotons(bool havePhotons);
	
	/** limit the step to a fraction of the mean free path
	 @param limit	fraction of the mean free path, should be between 0 and 1
	 */
	void setLimit(double limit);
	
	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);
	
	/** Apply a surface that confine the position dependent photon field
	 * @param surface closed surface to confine the grid nodes to be uploaded
	 */
	void setSurface(ref_ptr<Surface> surface);
	ref_ptr<Surface> getSurface() const;
	
	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	/** set a custom interaction rate
	 * With this function you can change the type of interaction rate,
	 * if you would for example like to change from a homogeneous to a position
	 * dependent interaction rate.
	 * @param intRates ref_ptr to a InteractionRates class
	 */
	void setInteractionRates(ref_ptr<InteractionRates> intRates);
	ref_ptr<InteractionRates> getInteractionRates() const;
	
	/** Loads the interaction rate
	 * This function loads the interaction rate from a given file/folder.
	 * @param path The name of the file/folder containing the interaction rates
	 */
	void initRate(std::string path);
	
	/** Loads the cumulative interaction rate
	 * This function loads the interaction rate from a given file/folder.
	 * @param path The name of the file/folder containing the interaction rates
	 */
	void initCumulativeRate(std::string path);


	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMINVERSECOMPTONSCATTERING_H

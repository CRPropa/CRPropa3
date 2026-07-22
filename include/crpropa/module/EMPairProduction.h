#ifndef CRPROPA_EMPAIRPRODUCTION_H
#define CRPROPA_EMPAIRPRODUCTION_H

#include <string>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Geometry.h"
#include "crpropa/InteractionRates.h"


namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */
 
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
 The surface is defined to include the nodes of the grid contained within.
 */
class EMPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField; 	// target photon field
	bool haveElectrons;					// add secondary electrons to simulation
	double limit;						// limit the step to a fraction of the mean free path
	double thinning;					// factor of the thinning (0: no thinning, 1: maximum thinning)
	ref_ptr<Surface> surface; // surface that includes the nodes in the photonField grid to be included
	std::string interactionTag = "EMPP";
	ref_ptr<InteractionRates> interactionRates;
public:
	/** Constructor
	 The object used to load, store and access to the interaction rates of the process is the interactionRates pointer.
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
	 @param surface  		Grid will be confined to `surface->distance(pos)<0`, so for example inside a closed surface
	 */
	EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1, ref_ptr<Surface> surface = nullptr);

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
	 * The rates are initialized only for distances smaller then 0, so `surface->distance(pos)<0`
	 * @param surface Grid will be confined to `surface->distance(pos)<0`, so for example inside a closed surface
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

	/**
	 * Get the interaction rate for a given energy, position, and redshift.
	 * For now, this function uses a redshift scaling factor for the interaction rate.
	 * Future releases will include a more accurate treatment of the redshift evolution of the photon field.
	 * @param E Energy of the primary particle
	 * @param position Position of the primary particle
	 * @param z Redshift of the primary particle
	 */
	double getRate(double E, const Vector3d &position = Vector3d(0.), double z = 0) const;
		
	void performInteraction(ref_ptr<Candidate> candidate) const;
	void process(ref_ptr<Candidate> candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_EMPAIRPRODUCTION_H

#ifndef CRPROPA_ELECTRONPAIRPRODUCTION_H
#define CRPROPA_ELECTRONPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class ElectronPairProduction
 @brief Electron-pair production of charged nuclei with background photons.

 This module simulates electron-pair production (known in the literature as
 the Bethe-Heitler process), in which a charged nucleus interacts with a
 background photon and part of its energy is converted into an
 electron-positron pair. Unlike photo-pion production, the nucleus itself
 is never reassigned a new particle ID -- only its energy is reduced -- and
 the fractional energy loss per interaction is small, so this module treats
 the process as a continuous, deterministic energy loss applied every step,
 rather than a discrete, randomly-sampled interaction.

 Energy loss rates are precomputed as a function of nucleus Lorentz factor
 for a given background photon field and loaded from a tabulated data file;
 above the tabulated range the rate is extrapolated using a power law,
 while below it no interaction is assumed. The tabulated rate (computed for
 protons) is rescaled for general nuclei using the nuclear charge and mass
 number, since the interaction depends on the square of the nuclear charge
 relative to its mass; the rate is further scaled with redshift to account
 for the evolving photon background. By default (limit = 0.1), the module
 limits each propagation step to 10% of the local energy-loss length,
 keeping the continuous-loss approximation accurate as conditions change
 along the trajectory.

 When secondary electron/positron pairs are activated, individual pairs are
 stochastically sampled from a tabulated electron energy spectrum (binned
 by nucleus Lorentz factor) until the sampled pair energies account for the
 total energy lost in that step, with the final pair accepted
 probabilistically if it would otherwise overshoot the remaining energy
 budget -- ensuring the discrete secondaries remain consistent, on average,
 with the continuous energy loss applied to the primary.

 @see PhotoPionProduction, NuclearDecay, Redshift
 */
class ElectronPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	std::vector<double> tabLossRate; /*< tabulated energy loss rate in [J/m] for protons at z = 0 */
	std::vector<double> tabLorentzFactor; /*< tabulated Lorentz factor */
	std::vector<std::vector<double> > tabSpectrum; /*< electron/positron cdf(Ee|log10(gamma)) for log10(Ee/eV)=7-24 in 170 steps and log10(gamma)=6-13 in 70 steps and*/
	double limit; ///< fraction of energy loss length to limit the next step
	bool haveElectrons; /*< if true, secondary electrons will be added to the simulation */
	std::string interactionTag = "EPP";

public:
	/**
	 * @brief Constructor for the Electron Pair Production
	 * 
	 * @param photonField 	target photon field
	 * @param haveElectrons If true, secondary electrons will be added to the simulation
	 * @param limit 		step size limit as fraction of mean free path
	 */
	ElectronPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons =
			false, double limit = 0.1);

	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary electrons are added to the simulation
	void setHaveElectrons(bool haveElectrons);
	
	/** Limit the propagation step to a fraction of the mean free path
	 * @param limit fraction of the mean free path
	 */
	void setLimit(double limit);
	
	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;
	
	void initRate(std::string filename);
	void initSpectrum(std::string filename);
	void process(Candidate *candidate) const;

	/**
	 Calculates the energy loss length 1/beta = -E dx/dE in [m]
	 @param	id		PDG particle ID
	 @param lf		Lorentz factor
	 @param z		redshift

	 The energy loss length is tabulated for protons against CMB and IRB.
	 Modification for nuclei and cosmological evolution of the photon background
	 is considered with (cf. 10.1016/j.astropartphys.2012.07.010, eq. 3 and 5)
	 beta_A,Z(E) = Z^2 / A * beta_p(E/A)
	 beta(E,z) = (1+z)^3 beta((1+z)E).
	 */
	double lossLength(int id, double lf, double z=0) const;
	
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_ELECTRONPAIRPRODUCTION_H

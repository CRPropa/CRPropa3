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

 This module simulates electron-pair production as a continuous energy loss.\n
 Several photon fields can be selected.\n
 The production of secondary e+/e- pairs and photons can by activated.\n
 By default, the module limits the step size to 10% of the energy loss length of the particle.
 */
class ElectronPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	std::vector<double> tabLossRate; /*< tabulated energy loss rate in [J/m] for protons at z = 0 */
	std::vector<double> tabLorentzFactor; /*< tabulated Lorentz factor */
	std::vector<std::vector<double> > tabSpectrum; /*< electron/positron cdf(Ee|log10(gamma)) for log10(Ee/eV)=7-24 in 170 steps and log10(gamma)=6-13 in 70 steps and*/
	double limit; ///< fraction of energy loss length to limit the next step
	bool haveElectrons;
	std::string interactionTag = "EPP";

public:
	ElectronPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons =
			false, double limit = 0.1);

	void setPhotonField(ref_ptr<PhotonField> photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);
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

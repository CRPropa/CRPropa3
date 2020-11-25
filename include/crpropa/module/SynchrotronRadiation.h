#ifndef CRPROPA_SYNCHROTRONRADIATION_H
#define CRPROPA_SYNCHROTRONRADIATION_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class SynchrotronRadiation
 @brief Synchrotron radiation of charged particles in magnetic fields.

 This module simulates the continuous energy loss of charged particles in magnetic fields, c.f. Jackson.
 The magnetic field is specified either by a MagneticField or by a RMS field strength value.
 The module limits the next step size to ensure a fractional energy loss dE/E < limit (default = 0.1).
 Optionally, synchrotron photons above a threshold (default E > 10^7 eV) are created as secondary particles.
 Note that the large number of secondary photons per propagation can cause memory problems.
 To mitigate this, use thinning. However, this still doesn't solve the problem completely.
 For this reason, a break-condition stops tracking secondary photons and reweights the current ones. 
 */
class SynchrotronRadiation: public Module {
private:
	ref_ptr<MagneticField> field; ///< MagneticField instance
	double Brms; ///< Brms value in case no MagneticField is specified
	double limit; ///< fraction of energy loss length to limit the next step
	double thinning; ///< thinning parameter for weighted-sampling (maximum 1, minimum 0)
	bool havePhotons; ///< flag for production of secondary photons
	int maximumSamples; ///< maximum number of samples of synchrotron photons (break condition; defaults to 100; 0 or <0 means no sampling)
	double secondaryThreshold; ///< threshold energy for secondary photons
	std::vector<double> tabx; ///< tabulated fraction E_photon/E_critical from 10^-6 to 10^2 in 801 log-spaced steps
	std::vector<double> tabCDF; ///< tabulated CDF of synchrotron spectrum


public:
	SynchrotronRadiation(ref_ptr<MagneticField> field, bool havePhotons = false, double thinning = 0, int nSamples = 0, double limit = 0.1);
	SynchrotronRadiation(double Brms = 0, bool havePhotons = false, double thinning = 0, int nSamples = 0, double limit = 0.1);
	void setField(ref_ptr<MagneticField> field);
	void setBrms(double Brms);	
	void setHavePhotons(bool havePhotons);
	void setThinning(double thinning);
	void setLimit(double limit);
	void setMaximumSamples(int nmax);
	void setSecondaryThreshold(double threshold);	
	ref_ptr<MagneticField> getField();
	double getBrms();
	bool getHavePhotons();
	double getThinning();
	double getLimit();
	int getMaximumSamples();
	double getSecondaryThreshold() const;
	void initSpectrum();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_SYNCHROTRONRADIATION_H

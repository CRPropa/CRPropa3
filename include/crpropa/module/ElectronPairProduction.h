#ifndef CRPROPA_ELECTRONPAIRPRODUCTION_H
#define CRPROPA_ELECTRONPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

namespace crpropa {

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
	PhotonField photonField;

	std::vector<double> tabLossRate; /*< tabulated energy loss rate in [J/m] for protons at z = 0 */
	std::vector<double> tabLorentzFactor; /*< tabulated Lorentz factor */

	std::vector<std::vector<double> > tabSpectrum; /*< tabulated electron spectrum dN/dE in [a.u.] */
	std::vector<double> tabE; /*< tabulated proton energy in [J], 70 steps from 10^15 - 10^22 eV */
	std::vector<double> tabEe; /*< edges of electron energy bins in [J], 171 steps from 10^6.95 - 10^23.95 eV */
	std::vector<double> tabEeWidth; /*< electron energy bin width in [J], 170 steps */
	double limit; ///< fraction of energy loss length to limit the next step
	bool haveElectrons;

public:
	ElectronPairProduction(PhotonField photonField = CMB, bool haveElectrons =
			false, double limit = 0.1);

	void setPhotonField(PhotonField photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initSpectrum(std::string filename);
	void process(Candidate *candidate) const;

	/**
	 Calculates the energy loss length 1/beta = -E dx/dE in [m]
	 @param	id		PDG particle ID
	 @param lf		Lorentz factor
	 @param z		redshift

	 The energy loss length is tabulated for protons against CMB, IRB and
	 CMB+IRB photon backgrounds.
	 Modification for nuclei and cosmological evolution of the photon background
	 is considered with (cf. 10.1016/j.astropartphys.2012.07.010, eq. 3 and 5)
	 beta_A,Z(E) = Z^2 / A * beta_p(E/A)
	 beta(E,z) = (1+z)^3 beta((1+z)E).
	 */
	double lossLength(int id, double lf, double z=0) const;
	void addElectrons(Candidate *candidate, double loss) const;
};

} // namespace crpropa

#endif // CRPROPA_ELECTRONPAIRPRODUCTION_H

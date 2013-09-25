#ifndef CRPROPA_ELECTRONPAIRPRODUCTION_H
#define CRPROPA_ELECTRONPAIRPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

namespace crpropa {

/**
 @class ElectronPairProduction
 @brief Electron-pair production of charged nuclei with background photons.

 This module simulates electron-pair production as a continuous energy loss.\n
 Several photon fields can be selected.
 Secondary electrons are not created.
 */
class ElectronPairProduction: public Module {
private:
	PhotonField photonField;
	std::vector<double> tabLossLength; /*< tabulated energy loss rate in [J/m] for protons at z = 0*/
	std::vector<double> tabLorentzFactor; /*< tabulated proton energy [J] */

public:
	ElectronPairProduction(PhotonField photonField = CMB_IRB);
	void setPhotonField(PhotonField photonField);
	void init();
	void init(std::string filename);
	void process(Candidate *candidate) const;

	/**
	 Calculates the inverse energy loss length beta = -1/E dE/dx in [1/m]
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
	double invLossLength(int id, double lf, double z) const;
};

} // namespace crpropa

#endif // CRPROPA_ELECTRONPAIRPRODUCTION_H

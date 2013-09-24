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
	std::vector<double> tabLossRate; /*< tabulated energy loss rate in [J/m] for protons at z = 0*/
	std::vector<double> tabEnergy; /*< tabulated proton energy [J] */

public:
	ElectronPairProduction(PhotonField photonField = CMB_IRB);
	void setPhotonField(PhotonField photonField);
	void init();
	void init(std::string filename);
	void process(Candidate *candidate) const;

	/**
	 Calculates the energy loss rate dE/dx in [J/m]
	 @param	id		PDG particle id
	 @param E		energy [J]
	 @param z		redshift

	 The energy loss rate b(E) = -dE/dt is tabulated for protons against
	 CMB, IRB and CMB+IRB photon backgrounds.

	 For nuclei this loss rate is modified as (cf. 10.1016/j.astropartphys.2012.07.010, eq. 5)
	 b_A,Z(E) = Z^2/A * b_p(E/A).

	 Cosmological evolution of the photon background is considered with (cf. 10.1103/PhysRevD.74.043005, eq. 5)
	 b(E,z) = (1+z)^2 b((1+z)E).
	 Note that the energy loss length beta(E) = -1/E dE/dt evolves as
	 beta(E,z) = (1+z)^3 beta((1+z)E).
	 */
	double lossRate(int id, double E, double z) const;
};

} // namespace crpropa

#endif // CRPROPA_ELECTRONPAIRPRODUCTION_H

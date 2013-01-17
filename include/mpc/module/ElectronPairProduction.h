#ifndef MPC_ELECTRONPAIRPRODUCTION_H_
#define MPC_ELECTRONPAIRPRODUCTION_H_

#include "mpc/Module.h"
#include "mpc/PhotonBackground.h"

namespace mpc {

/**
 @class ElectronPairProduction
 @brief Electron-pair production of charged nuclei with background photons.

 This module simulates electron-pair production as a continuous energy loss.\n
 Several photon fields can be selected. They are considered homogeneous and evolving as the CMB.
 */
class ElectronPairProduction: public Module {
private:
	PhotonField photonField;
	std::vector<double> lossRate; // energy loss rate in [J/m]
	std::vector<double> energy;

public:
	ElectronPairProduction(PhotonField photonField = CMB_IRB);
	void setPhotonField(PhotonField photonField);
	void init();
	void init(std::string filename);
	void process(Candidate *candidate) const;

	/**
	 Calculates the energy loss length 1/E dE/dx in [m]
	 @param	id		PDG particle id
	 @param energy	particle energy [J]
	 */
	double energyLossLength(int id, double energy);
};

} // namespace mpc

#endif // MPC_ELECTRONPAIRPRODUCTION_H_

#ifndef MPC_PHOTOPIONPRODUCTION_H_
#define MPC_PHOTOPIONPRODUCTION_H_

#include "mpc/module/StochasticInteraction.h"
#include "mpc/Random.h"

#include <vector>

namespace mpc {

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.

 This module simulates photo-hadronic interactions of nuclei with background photons.\n
 Several photon fields can be selected. They are considered as homogeneous and evolving as the CMB.\n
 */
class PhotoPionProduction: public StochasticInteraction {
protected:
	PhotonField photonField;
	std::vector<double> pRate; // interaction rate in [1/m] for protons
	std::vector<double> nRate; // interaction rate in [1/m] for neutrons
	std::vector<double> energy; // energy in [J]

public:
	PhotoPionProduction(PhotonField photonField = CMB);
	void setPhotonField(PhotonField photonField);
	void init();
	void init(std::string filename);
	bool setNextInteraction(Candidate *candidate,
			InteractionState &interaction) const;
	void performInteraction(Candidate *candidate) const;

	/**
	 Calculates the energy loss length 1/E dE/dx in [m]
	 @param	id		PDG particle id
	 @param energy	particle energy [J]
	 */
	double energyLossLength(int id, double energy);
};

/**
 @class SophiaPhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons using SOPHIA.

 This module simulates photo-hadronic interactions of nuclei with background photons.\n
 Several photon fields can be selected. They are considered as homogeneous and evolving as the CMB.\n
 The interaction itself is simulated with SOPHIA.\n
 Electromagnetic particles, neutrinos and antiparticles as secondaries from these interactions can be switched on independently.
 */
class SophiaPhotoPionProduction: public PhotoPionProduction {
private:
	bool havePhotons;
	bool haveNeutrinos;
	bool haveAntiNucleons;

public:
	SophiaPhotoPionProduction(PhotonField photonField = CMB, bool photons =
			false, bool neutrinos = false, bool antiNucleons = false);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void setHaveAntiNucleons(bool b);
	void performInteraction(Candidate *candidate) const;
};

} // namespace mpc

#endif /* MPC_PHOTOPIONPRODUCTION_H_ */

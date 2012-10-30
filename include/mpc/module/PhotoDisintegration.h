#ifndef PHOTODISINTEGRATION_H_
#define PHOTODISINTEGRATION_H_

#include "mpc/module/StochasticInteraction.h"
#include "mpc/Random.h"

#include <vector>

namespace mpc {

/**
 @class PhotoDisintegration
 @brief Photo-disintegration of nuclei with background photons.

 This module simulates photo-disintegration of nuclei with background photons.\n
 Background photon fields are considered as homogeneous and evolving as the CMB.\n
 */
class PhotoDisintegration: public StochasticInteraction {
private:
	PhotonField photonField;
	struct PDMode {
		int channel; // number of emitted (n, p, H2, H3, He3, He4)
		std::vector<double> rate; // disintegration rate [1/m]
	};
	std::vector<std::vector<PDMode> > pdTable; // pdTable[Z * 31 + N] = vector<PDmode>

public:
	PhotoDisintegration(PhotonField photonField = CMB);
	void init(PhotonField photonField);
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

} // namespace mpc

#endif /* PHOTODISINTEGRATION_H_ */

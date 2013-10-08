#ifndef CRPROPA_PHOTOPIONPRODUCTION_H
#define CRPROPA_PHOTOPIONPRODUCTION_H

#include "crpropa/module/StochasticInteraction.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.
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

	bool randomInteraction(Candidate *candidate,
			InteractionState &interaction) const;

	void performInteraction(Candidate *candidate,
			InteractionState &interaction) const;

	double nucleiModification(int A, int X) const;

	/**
	 Calculates the energy loss length 1/E dE/dx in [m]. This is not used in the simulation.
	 @param	id		PDG particle id
	 @param energy	particle energy [J]
	 */
	double energyLossLength(int id, double energy);
};

/**
 @class SophiaPhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons using SOPHIA.

 This module simulates photo-hadronic interactions of nuclei with background photons.\n
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
	void performInteraction(Candidate *candidate,
			InteractionState &interaction) const;
};

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTION_H

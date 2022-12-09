#ifndef CRPROPA_PHOTODISINTEGRATION_H
#define CRPROPA_PHOTODISINTEGRATION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>
#include <map>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */
 
/**
 @class PhotoDisintegration
 @brief Photodisintegration of nuclei by background photons.
 */
class PhotoDisintegration: public Module {
private:
	ref_ptr<PhotonField> photonField;
	double limit; // fraction of mean free path for limiting the next step
	bool havePhotons;
	std::string interactionTag = "PD";

	struct Branch {
		int channel; // number of emitted (n, p, H2, H3, He3, He4)
		std::vector<double> branchingRatio; // branching ratio as function of nucleus Lorentz factor
	};

	struct PhotonEmission {
		double energy; // energy of emitted photon [J]
		std::vector<double> emissionProbability; // emission probability as function of nucleus Lorentz factor
	};

	std::vector<std::vector<double> > pdRate; // pdRate[Z * 31 + N] = total interaction rate
	std::vector<std::vector<Branch> > pdBranch; // pdTable[Z * 31 + N] = branching ratios
	mutable std::map<int, std::vector<PhotonEmission> > pdPhoton; // map of emitted photon energies and photon emission probabilities

	static const double lgmin; // minimum log10(Lorentz-factor)
	static const double lgmax; // maximum log10(Lorentz-factor)
	static const size_t nlg; // number of Lorentz-factor steps

public:
	/** Constructor.
	 @param photonField		target photon field
	 @param havePhotons		if true, add secondary photons as candidates
	 @param limit			step size limit as fraction of mean free path
	 */
	PhotoDisintegration(ref_ptr<PhotonField> photonField, bool havePhotons = false, double limit = 0.1);

	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary photons are added to the simulation
	void setHavePhotons(bool havePhotons);

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
	void initBranching(std::string filename);
	void initPhotonEmission(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;

	/**
	 Calculates the loss length E dx/dE in [m] physical distance.
	 This is not used in the simulation.
	 @param	id		PDG particle id
	 @param gamma	Lorentz factor of particle
	 @param z		redshift
	 @returns E dx/dE [in meters]
	 */
	double lossLength(int id, double gamma, double z = 0);
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PHOTODISINTEGRATION_H

#ifndef CRPROPA_NUCLEARDECAY_H
#define CRPROPA_NUCLEARDECAY_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class NuclearDecay
 @brief Nuclear decay of unstable nuclei.

 This module simulates the nuclear decay of unstable nuclei using data from NuDat2.
 All decay modes are considered: alpha, beta+- and gamma decay, as well as proton- and neutron dripping.
 The resulting non-hadronic secondary particles (e+, e-, neutrinos, gamma) can optionally be created.

 For details on the preprocessing of the NuDat2 data refer to "CRPropa3-data/calc_decay.py".
 */
class NuclearDecay: public Module {
private:
	double limit;
	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;
	struct DecayMode {
		int channel; // (#beta- #beta+ #alpha #proton #neutron)
		double rate; // decay rate in [1/m]
		std::vector<double> energy; // photon energies of ensuing gamma decays
		std::vector<double> intensity; // probabilities of ensuing gamma decays
	};
	std::vector<std::vector<DecayMode> > decayTable; // decayTable[Z * 31 + N] = vector<DecayMode>
	std::string interactionTag = "ND";

public:
	/** Constructor.
	 @param electrons		if true, add secondary photons as candidates
	 @param photons			if true, add secondary photons as candidates
	 @param neutrinos		if true, add secondary neutrinos as candidates
	 @param limit			step size limit as fraction of mean free path
	 */
	NuclearDecay(bool electrons = false, bool photons = false, bool neutrinos = false, double limit = 0.1);

	/** Limit the propagation step to a fraction of the mean free path
	 * @param limit fraction of the mean free path
	 */
	void setLimit(double limit);

	// decide if secondary electrons are added to the simulation	
	void setHaveElectrons(bool b);

	// decide if secondary photons are added to the simulation	
	void setHavePhotons(bool b);

	// decide if secondary neutrinos are added to the simulation	
	void setHaveNeutrinos(bool b);

	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	void gammaEmission(Candidate *candidate, int channel) const;
	void betaDecay(Candidate *candidate, bool isBetaPlus) const;
	void nucleonEmission(Candidate *candidate, int dA, int dZ) const;

	/**
	 Return the mean free path.
	 This is not used in the simulation.
	 @param id      PDG particle id
	 @param gamma   Lorentz factor of particle
	 @returns The mean free path [in meters]
	 */
	double meanFreePath(int id, double gamma);
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_NUCLEARDECAY_H

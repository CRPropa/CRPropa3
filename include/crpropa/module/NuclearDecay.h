#ifndef CRPROPA_NUCLEARDECAY_H
#define CRPROPA_NUCLEARDECAY_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {

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

public:
	NuclearDecay(bool electrons = false, bool photons = false, bool neutrinos = false, double limit = 0.1);
	void setLimit(double limit);
	void setHaveElectrons(bool b);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	void gammaEmission(Candidate *candidate, int channel) const;
	void betaDecay(Candidate *candidate, bool isBetaPlus) const;
	void nucleonEmission(Candidate *candidate, int dA, int dZ) const;

	/**
	 Return the mean free path in [m].
	 This is not used in the simulation.
	 @param id      PDG particle id
	 @param gamma   Lorentz factor of particle
	 */
	double meanFreePath(int id, double gamma);
};

} // namespace crpropa

#endif // CRPROPA_NUCLEARDECAY_H

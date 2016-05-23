#ifndef CRPROPA_NUCLEARDECAY_H
#define CRPROPA_NUCLEARDECAY_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {

/**
 @class NuclearDecay
 @brief Nuclear decay of unstable nuclei.

 This module simulates the nuclear decay of unstable nuclei using data from NuDat2.
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
		std::vector<double> energy;
		std::vector<double> intensity;
	};
	std::vector<std::vector<DecayMode> > decayTable; // decayTable[Z * 31 + N] = vector<DecayMode>

public:
	NuclearDecay(bool electrons = false, bool photons = false, bool neutrinos = false, double limit =
			0.1);
	void setLimit(double limit);
	void setHaveElectrons(bool b);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	void gammaEmission(Candidate *candidate, int channel, double &Egamma) const;
	void betaDecay(Candidate *candidate, bool isBetaPlus, double Egamma) const;
	void nucleonEmission(Candidate *candidate, int dA, int dZ) const;

    /**
     Calculates the loss length E dx/dE in [m].
     This is not used in the simulation.
     @param id      PDG particle id
     @param gamma   Lorentz factor of particle
     @param z       redshift
     */
    double lossLength(int id, double gamma, double z = 0);
};

} // namespace crpropa

#endif // CRPROPA_NUCLEARDECAY_H

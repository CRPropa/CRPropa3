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
	bool haveNeutrinos;
	struct DecayMode {
		int channel; // (#beta- #beta+ #alpha #proton #neutron)
		double rate; // decay rate in [1/m]
	};
	std::vector<std::vector<DecayMode> > decayTable; // decayTable[Z * 31 + N] = vector<DecayMode>
	std::vector<double> tBeta; // electron kinetic energy [J] in neutron decays
	std::vector<double> cdfBeta; // cumulative distribution function for the electron kinetic energy [J] in neutron decays

public:
	NuclearDecay(bool electrons = false, bool neutrinos = false, double limit =
			0.1);
	void setLimit(double limit);
	void setHaveElectrons(bool b);
	void setHaveNeutrinos(bool b);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	void betaDecay(Candidate *candidate, bool isBetaPlus) const;
	void nucleonEmission(Candidate *candidate, int dA, int dZ) const;
};

} // namespace crpropa

#endif // CRPROPA_NUCLEARDECAY_H

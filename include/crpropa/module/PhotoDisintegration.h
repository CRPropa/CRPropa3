#ifndef CRPROPA_PHOTODISINTEGRATION_H
#define CRPROPA_PHOTODISINTEGRATION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {

/**
 @class PhotoDisintegration
 @brief Photo-disintegration of nuclei with background photons.
 */
class PhotoDisintegration: public Module {
private:
	PhotonField photonField;
	double limit; // fraction of mean free path for limiting the next step

	struct Branch {
		int channel; // number of emitted (n, p, H2, H3, He3, He4)
		std::vector<double> branchingRatio; // branching ratio as function of nucleus Lorentz factor
	};
	std::vector<std::vector<Branch> > pdBranch; // pdTable[Z * 31 + N] = vector<Branch>
	std::vector<std::vector<double> > pdRate; // pdRate[Z * 31 + N] = total disintegration rate

	static const double lgmin; // minimum log10(Lorentz-factor)
	static const double lgmax; // maximum log10(Lorentz-factor)
	static const size_t nlg; // number of Lorentz-factor steps

public:
	PhotoDisintegration(PhotonField photonField = CMB, double limit = 0.1);

	void setPhotonField(PhotonField photonField);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initBranching(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;

	/**
	 Calculates the loss length E dx/dE in [m] physical distance.
	 This is not used in the simulation.
	 @param	id		PDG particle id
	 @param gamma	Lorentz factor of particle
	 @param z		redshift
	 */
	double lossLength(int id, double gamma, double z = 0);
};

} // namespace crpropa

#endif // CRPROPA_PHOTODISINTEGRATION_H

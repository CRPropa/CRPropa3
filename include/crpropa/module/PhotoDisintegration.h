#ifndef CRPROPA_PHOTODISINTEGRATION_H
#define CRPROPA_PHOTODISINTEGRATION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>
#include <map>

namespace crpropa {

/**
 @class PhotoDisintegration
 @brief Photo-disintegration of nuclei with background photons.
 */
class PhotoDisintegration: public Module {
private:
	PhotonField photonField;
	double limit; // fraction of mean free path for limiting the next step
	bool havePhotons;

	struct Branch {
		int channel; // number of emitted (n, p, H2, H3, He3, He4)
		std::vector<double> branchingRatio; // branching ratio as function of nucleus Lorentz factor
	};
	struct PhotonEmission {
		double energy; // energy of emitted photon [J]
		std::vector<double> emissionProbability; // emission probability as function of nucleus Lorentz factor
	};
	std::vector<std::vector<Branch> > pdBranch; // pdTable[Z * 31 + N] = vector<Branch>
	std::vector<std::vector<double> > pdRate; // pdRate[Z * 31 + N] = total disintegration rate
	mutable std::map<int, std::vector<PhotonEmission> > pdPhoton; // map of emitted photon energies, photon emission probability as function of gamma from 1e6 to 1e14 within 201 logspaced steps for each combination of mother and daughter isotopes 
	std::vector<std::vector<double> > elasticRate; // elatsicRate[Z * 31 + N] = elastic scattering rate
//	std::vector<std::vector<double> > elasticCDF; //CDF as function of background photon energy in nucleus restframe from 1e3 eV to 2.63e8 eV in 513 steps. For each gamma factor from 1e6 to 1e14 within 201 steps a average CDF is calculated over all isotopes. TODO: rewrite title
	std::vector<std::vector<double> > elasticCDF; //CDF as function of background photon energy in nucleus restframe from 2e3 eV to 2.63e8 eV in 513 steps. For each gamma factor from 1e6 to 1e14 within 201 steps a average CDF is calculated over all isotopes.
	std::vector<double> lgElastic; // log10 gamma factor for tabulated CDFs

	static const double lgmin; // minimum log10(Lorentz-factor)
	static const double lgmax; // maximum log10(Lorentz-factor)
	static const size_t nlg; // number of Lorentz-factor steps
	static const double epsMin; // minimum log10(eps / J)
	static const double epsMax; // maximum log10(eps / J)
	static const size_t neps; // number of eps steps
	std::vector<double> tabEps; // background photon energies [J] for elastic scattering in nucleus rest frame from 2 keV to 264 MeV in 513 steps

public:
	PhotoDisintegration(PhotonField photonField = CMB, bool havePhotons = false, double limit = 0.1);

	void setPhotonField(PhotonField photonField);
	void setHavePhotons(bool havePhotons);
	void setLimit(double limit);

	void initRate(std::string filename);
	void initBranching(std::string filename);
	void initPhotonEmission(std::string filename);
	void initElastic(std::string filename);
	void initElasticCDF(std::string filename);

	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	void performElasticScattering(Candidate *candidate, std::vector<double> cdf) const;

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

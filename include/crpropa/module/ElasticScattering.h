#ifndef CRPROPA_ELASTICSCATTERING_H
#define CRPROPA_ELASTICSCATTERING_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {

/**
 @class ElasticScattering
 @brief Elastic scattering of background photons on cosmic-ray nuclei.
 */
class ElasticScattering: public Module {
private:
	ref_ptr<PhotonField> photonField;

	std::vector<double> tabRate; // elastic scattering rate
	std::vector<std::vector<double> > tabCDF; // CDF as function of background photon energy
	std::string interactionTag = "ES";

	static const double lgmin; // minimum log10(Lorentz-factor)
	static const double lgmax; // maximum log10(Lorentz-factor)
	static const size_t nlg;   // number of Lorentz-factor steps
	static const double epsmin; // minimum log10(eps / J)
	static const double epsmax; // maximum log10(eps / J)
	static const size_t neps;   // number of eps steps

public:
	/** Constructor
	 @param photonField		target photon field
	 */
	ElasticScattering(ref_ptr<PhotonField> photonField);
	void initRate(std::string filename);
	void initCDF(std::string filename);
	void setPhotonField(ref_ptr<PhotonField> photonField);
	void process(Candidate *candidate) const;
	
	std::string getInteractionTag() const;
	void setInteractionTag(std::string tag);
};

} // namespace crpropa

#endif // CRPROPA_ELASTICSCATTERING_H

#ifndef CRPROPA_ELASTICSCATTERING_H
#define CRPROPA_ELASTICSCATTERING_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {

/**
 @class ElasticScattering
 @brief Elastic scattering of background photons on cosmic ray nuclei.
 */
class ElasticScattering: public Module {
private:
    PhotonField photonField;
    ScalarGrid4d spaceTimeGrid;
    ScalarGrid spaceGrid;
    std::string tag;

    std::vector<double> tabRate; // elastic scattering rate
    std::vector<std::vector<double> > tabCDF; // CDF as function of background photon energy

    static const double lgmin; // minimum log10(Lorentz-factor)
    static const double lgmax; // maximum log10(Lorentz-factor)
    static const size_t nlg;   // number of Lorentz-factor steps
    static const double epsmin; // minimum log10(eps / J)
    static const double epsmax; // maximum log10(eps / J)
    static const size_t neps;   // number of eps steps

public:
    ElasticScattering(PhotonField photonField,
                      std::string tag = "ElScatt"
    );

    ElasticScattering(PhotonField photonField,
                      ScalarGrid4d spaceTimeGrid,
                      std::string tag = "ElScatt"
    );

    ElasticScattering(PhotonField photonField,
                      ScalarGrid spaceGrid,
                      std::string tag = "ElScatt"
    );

    void initRate(std::string filename);
    void initCDF(std::string filename);
    void setPhotonField(PhotonField photonField);
    void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_ELASTICSCATTERING_H

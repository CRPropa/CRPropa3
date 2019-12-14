#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include <string>

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */
// Photon fields
// The default IRB model is that of Kneiske et al. 2004
enum PhotonField {
	CMB,
	IRB,  // same as IRB_Kneiske04
	IRB_Kneiske04,
	IRB_Stecker05,
	IRB_Franceschini08,
	IRB_Finke10,
	IRB_Dominguez11,
	IRB_Gilmore12,
	IRB_Stecker16_upper,
	IRB_Stecker16_lower,
	URB_Protheroe96
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

// Returns a string representation of the field
std::string photonFieldName(PhotonField photonField);


/*
	methods related to photon sampling as done in SOPHIA.
	only needed in PhotoPionProduction for this implementation.
*/

class Photon_Field {
 public:
 	Photon_Field();
    explicit Photon_Field(int bgFlag);
    double sample_eps(bool onProton, double E_in, double z_in) const;
 protected:
    int bgFlag;
    double getPhotonDensity(double eps, double z_in) const;
    double gaussInt(std::string type, double lowerLimit, double upperLimit, bool onProton, double E_in, double z_in) const;
    double functs(double s, bool onProton) const;
    double prob_eps(double eps, bool onProton, double E_in, double z_in) const;
    double crossection(double eps, bool onProton) const;
        double Pl(double x, double xth, double xmax, double alpha) const;
        double Ef(double x, double th, double w) const;
        double breitwigner(double sigma_0, double Gamma, double DMM, double epsPrime, bool onProton) const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H

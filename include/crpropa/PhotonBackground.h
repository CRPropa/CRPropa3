#ifndef CRPROPA_PHOTONBACKGROUND_H
#define CRPROPA_PHOTONBACKGROUND_H

#include <crpropa/Vector3.h>
#include <crpropa/Grid.h>

#include <iostream>  // cout, endl
#include <cmath>  // sqrt, pow
#include <string>
#include <fstream>  // write to file
#include <algorithm>  // max_element
#include <limits>  // for ::max
#include <vector>

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */
// Photon fields
// The default IRB model is that of Kneiske et al. 2004
// The slots PF1 to PF8 may be used for custom photon fields
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
	URB_Protheroe96,
	PF1, PF2, PF3, PF4,  // customizable
	PF5, PF6, PF7, PF8,  // field slots
};

// Returns overall comoving scaling factor
double photonFieldScaling(PhotonField photonField, double z);

// Returns a string representation of the field
std::string photonFieldName(PhotonField photonField);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
// custom photon field methods
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Photon_Field {
 public:
    explicit Photon_Field(std::string fieldPath);
    Photon_Field();
    double sample_eps(bool onProton, double E_in, double z_in) const;

 private:
    void init(std::string fieldPath);
        std::vector<double> energy;
        std::vector< std::vector<double> > dn_deps;
        std::vector<double> redshift;
    double get_photonDensity(double eps, int z_pos) const;
    double gaussInt(std::string type, double lowerLimit, double upperLimit, bool onProton, double E_in, int z_pos) const;
    double functs(double s, bool onProton) const;
    double prob_eps(double eps, bool onProton, double E_in, int z_pos) const;
    double crossection(double eps, bool onProton) const;
        double Pl(double, double, double, double) const;
        double Ef(double, double, double) const;
        double breitwigner(double, double, double, double, bool onProton) const;
        double singleback(double) const;
        double twoback(double) const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PHOTONBACKGROUND_H

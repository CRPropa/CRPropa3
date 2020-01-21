#ifndef CRPROPA_TD13FIELD_H
#define CRPROPA_TD13FIELD_H

#include <crpropa/magneticField/MagneticField.h>
#include <crpropa/Vector3.h>

#include <vector>

using namespace std;
using namespace crpropa;

/* 
 @class TD13Field
 @brief TD13Field turbulent magnetic field 

 Implementation of the turbulent magnetic field from Tautz and Dosch 2013
 "On numerical turbulence generation for test-particle simulations"
 in Physics of Plasmas
 doi: 10.1063/1.4789861
 see https://ui.adsabs.harvard.edu/abs/2013PhPl...20b2302T/abstract
 */

class TD13Field : public MagneticField {
private:
	double spec_Lmax, spec_Lmin; 
    double spec_q, spec_s; 
    int Nmodes; 

    vector<double> k_n;
    vector<double> Ak_n;
    vector<double> eta_n; // = cos(theta_n)
    vector<double> sqrt_eta_n; // = sqrt(1-eta**2)
    vector<double> cos_phi_n;
    vector<double> sin_phi_n;
    vector<double> phase_n;
    vector<Vector3d> Xi_n;

public:
    /**Constructor
       @param B_0           Magnetic field strength at reference level
       @param Lmin, Lmax    Magnetic field minimum and maximum wave length
       @param s, q          Magnetic field spectral indexes : default values (s=5/3, q=0) for a Kolmogorov spectrum 
       @param Nm            Number of Fourier modes
       @param seed          Seed for the random number generator, if 0, not set
     */
    TD13Field(double B_0, double Lmin, double Lmax, double s=5./3., double q=0, int Nm=64, int seed=42);

    Vector3d getField(const Vector3d &pos) const;

    /**@brief       compute the magnetic field coherence length according to the formula in  Harari et Al JHEP03(2002)045  
     * @return Lc   coherence length of the magnetic field
     */
    double getLc() const;
};

#endif // CRPROPA_TD13FIELD_H

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

=======
#include <vector>
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {

std::vector<double> logspace(double start, double stop, size_t N) {

  double delta = stop - start;
  std::vector<double> values = std::vector<double>(N, 0.);
  for (int i=0; i<N; i++) {
    values[i] = pow(10, ((double) i) / ((double) (N-1)) * delta + start);
  }
  return values;
} 

/**
 @class TD13Field
 @brief Interpolation-free turbulent magnetic field based on the TD13 paper

 blah blah blah
 */
class TD13Field: public MagneticField {
private:

  std::vector<Vector3d> xi; //TODO: this is actually psi (as defined in the paper), because I'm stupid. I don't think that's a problem, but I should probably still change it.
  std::vector<Vector3d> kappa;
  std::vector<double> phi;
  std::vector<double> costheta;
  std::vector<double> beta;
  std::vector<double> k;
  std::vector<double> Ak;

  double gamma;
  double Nm;

public:
  /** Constructor
      @param kmin wave number of the mode with the largest wavelength to be included in the spectrum
      @param kmax wave number of the mode with the smallest wavelength to be included in the spectrum
      @param gamma spectral index
      @param Nm number of wavemodes that will be used when computing the field. A higher value will give a more accurate representation of the turbulence, but increase the runtime for getField.
*/
  TD13Field(double kmin, double kmax, double gamma, double Nm);

  /**
     Theoretical runtime is O(Nm).
*/
  Vector3d getField(const Vector3d& pos) const;
};

} // namespace crpropa

#endif // CRPROPA_TD13FIELD_H

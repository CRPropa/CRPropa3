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

#include <immintrin.h>

namespace crpropa {

std::vector<double> logspace(double start, double stop, size_t N);

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

  double gamma;
  double Nm;


  int avx_Nm;
  int align_offset;

  static const int ixi0 = 0;
  static const int ixi1 = 1;
  static const int ixi2 = 2;
  static const int ikappa0 = 3;
  static const int ikappa1 = 4;
  static const int ikappa2 = 5;
  static const int iAk = 6;
  static const int ik = 7;
  static const int ibeta = 8;
  static const int itotal = 9;

public:
  //DEBUG put these in public so I can read them in python
  std::vector<float> avx_data;
  std::vector<double> Ak;
  std::vector<double> k;
  /** Constructor
      @param root mean square field strength for generated field
      @param kmin wave number of the mode with the largest wavelength to be included in the spectrum
      @param kmax wave number of the mode with the smallest wavelength to be included in the spectrum
      @param gamma spectral index
      @param bendoverScale the turbulence bend-over scale, used to scale the wavenumbers. As per the TD13 paper, this is set to 0.03AU by default.
      @param Nm number of wavemodes that will be used when computing the field. A higher value will give a more accurate representation of the turbulence, but increase the runtime for getField.
      @param seed can be used to seed the random number generator used to generate the field. This works just like in initTurbulence: a seed of 0 will lead to a randomly initialized RNG.
*/
  TD13Field(double Brms, double kmin, double kmax, double gamma,
	    double bendoverScale = 4.5e9, int Nm = 100, int seed = 0);

  // TODO: bendoverScale: figure out how to improve this, b/c it takes up space in the constructor arguments
  // (Lukas's suggestion was to use a setter for this, but it turns out that's not that easy
  // to do, b/c the bendoverScale is _only_ used in the constructor. We could only use the setter if the field
  // generation was moved out of the constructor, but I'm not sure that's desirable.)

  /**
     Theoretical runtime is O(Nm).
*/
  Vector3d getField(const Vector3d& pos) const;

  // versions:
  // 4: introduce field versioning;
  //    add volume correction factor   
  static int fieldVersion() {
    return 4;
  }
};

} // namespace crpropa

#endif // CRPROPA_TD13FIELD_H

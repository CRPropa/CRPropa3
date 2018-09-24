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

std::vector<double> logspace(double start, double stop, size_t N);

/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class TD13Field
 @brief Interpolation-free turbulent magnetic field based on the GJ99 and TD13 papers

 ## Overview
 This class provides a turbulent magnetic field that is generated as described by [(Giacalone and Jokipii, 1999)][gj99] and [(Tautz and Dosch, 2013)][td13]. Instead of using an inverse Fourier transform to generate the field on a grid -- which then needs to be interpolated to obtain in-between values -- this method only generates the wave modes making up the turbulent magnetic field ahead of time. At run time, when the field's value at a particular position is required, these plane waves are then evaluated analytically at that position. This guarantees that the resulting field is completely free of divergence, reproduces the mean field strength accurately, and does not suffer from other interpolation-induced problems. The disadvantage is that the number of wave modes is drastically smaller when compared with initTurbulence, which might have physical ramifications on the particles propagating through the field. Additionally, the implementation is somewhat slower.

 ## Using the SIMD optimization
 In order to mitigate some of the performance impact that is inherent in this method of field generation, an optimized version utilizing data-level parallelism through SIMD instructions is provided. More specifically, this implementation uses the x86 extensions SSE1, SSE2 and SSE3, with SSE4.1 being optional. In order to use this optimized version, three conditions need to be met:

1. The vectorized math library [SLEEF](sleef.org) needs to be present in CRPropa's prefix path. This library provides a fast, vectorized cosine function required by getField(). When running CMake, the output will indicate whether SLEEF was found or not.
2. The `USE_SIMD` option needs to be explicitly enabled in CMake. Currently, this sets GCC flags that tell the compiler to allow SIMD instructions.
3. Finally, the CPU that will actually run the code needs to support the abovementioned extensions: SSE1 through SSE3. These extensions are relatively old and quite common, so this should not be a problem in practice.

 ### Installation of SLEEF
 To install SLEEF, clone its [GitHub Repository](https://github.com/shibatch/sleef/) and build the library. SLEEF uses the same build system as CRPropa (CMake), so you should be familiar with the required steps: Create a build directory, `cd` into it, and run cmake (you can also refer to [SLEEF's installation instructions](https://sleef.org/compile.xhtml)). Make sure to set `-DCMAKE_INSTALL_PREFIX` appropriately; if you're following [CRPropa's virtualenv installation guide](https://github.com/CRPropa/CRPropa3/wiki/Installation#install-from-source-in-python-virtualenv---a-detailed-guide-recommended), set it to the same value as for the CRPropa cmake call (i.e. `-DCMAKE_INSTALL_PREFIX=$CRPROPA_DIR`).

[gj99]: https://doi.org/10.1086/307452
[td13]: https://doi.org/10.1063/1.4789861
 */
class TD13Field: public MagneticField {
private:

  std::vector<Vector3d> xi; //TODO: this is actually psi (as defined in the paper), because I'm stupid. I don't think that's a problem, but I should probably still change it.
  std::vector<Vector3d> kappa;
  std::vector<double> phi;
  std::vector<double> costheta;
  std::vector<double> beta;

  std::vector<float> avx_data;
  std::vector<double> Ak;
  std::vector<double> k;

  int Nm;

  int avx_Nm;
  int align_offset;

  // iAxi is a combined array containing the product of A * xi
  static const int iAxi0 = 0;
  static const int iAxi1 = 1;
  static const int iAxi2 = 2;
  // ikkappa is a combined array containing the product of k * kappa
  static const int ikkappa0 = 3;
  static const int ikkappa1 = 4;
  static const int ikkappa2 = 5;
  static const int ibeta = 6;
  static const int itotal = 7;

public:
  /**
      Create a new instance of TD13Field with the specified parameters. This generates all of the wavemodes according to the given parameters.
      @param Brms root mean square field strength for generated field
      @param kmin wave number of the mode with the largest wavelength to be included in the spectrum
      @param kmax wave number of the mode with the smallest wavelength to be included in the spectrum
      @param gamma the spectral index, as defined in the TD13 paper. Usually, you'd want to use 5/3 here, which will be comparable to passing -11/3 to `initTurbulence`.
      @param bendoverScale the turbulence bend-over scale, used to scale the wavenumbers. As per the TD13 paper, this is set to 0.03AU by default.
      @param Nm number of wavemodes that will be used when computing the field. A higher value will give a more accurate representation of the turbulence, but increase the runtime for getField.
      @param seed can be used to seed the random number generator used to generate the field. This works just like in initTurbulence: a seed of 0 will lead to a randomly initialized RNG.
*/
  TD13Field(double Brms, double kmin, double kmax, double s = 5/3.,
	    double bendoverScale = 4.5e9, int Nm = 100, int seed = 0);

  // TODO: bendoverScale: figure out how to improve this, b/c it takes up space in the constructor arguments
  // (Lukas's suggestion was to use a setter for this, but it turns out that's not that easy
  // to do, b/c the bendoverScale is _only_ used in the constructor. We could only use the setter if the field
  // generation was moved out of the constructor, but I'm not sure that's desirable.)

  /**
     Evaluates the field at the given position.

     Theoretical runtime is O(Nm), where Nm is the number of wavemodes.
*/
  Vector3d getField(const Vector3d& pos) const;

  // versions:
  // 4: introduce field versioning;
  //    add volume correction factor   
  static int fieldVersion() {
    return 4;
  }
};

/** @} */

} // namespace crpropa

#endif // CRPROPA_TD13FIELD_H

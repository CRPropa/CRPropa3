#ifndef CRPROPA_PLANEWAVETURBULENCE_H
#define CRPROPA_PLANEWAVETURBULENCE_H

#include "crpropa/Grid.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include <vector>

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class PlaneWaveTurbulence
 @brief Interpolation-free turbulent magnetic field based on the GJ99 and TD13
papers

 ## Overview
 This class provides a turbulent magnetic field that is generated as described
by [(Giacalone and Jokipii, 1999)][gj99] and [(Tautz and Dosch, 2013)][td13].
Instead of using an inverse Fourier transform to generate the field on a grid --
which then needs to be interpolated to obtain in-between values -- this method
only generates the wave modes making up the turbulent magnetic field ahead of
time. At run time, when the field's value at a particular position is required,
these plane waves are then evaluated analytically at that position. This
guarantees that the resulting field is completely free of divergence, reproduces
the mean field strength accurately, and does not suffer from other
interpolation-induced problems. The disadvantage is that the number of wave
modes is drastically smaller when compared with initTurbulence, which might have
physical ramifications on the particles propagating through the field.
Additionally, the implementation is somewhat slower.

 ## Using the SIMD optimization
 In order to mitigate some of the performance impact that is inherent in this
method of field generation, an optimized version utilizing data-level
parallelism through SIMD instructions is provided. More specifically, this
implementation uses the AVX extension to the x86 instruction set. In order to
use this optimized version, three conditions need to be met:

1. The `USE_SIMD` option needs to be explicitly enabled in CMake. Currently,
this sets GCC flags that tell the compiler to allow SIMD instructions.
2. You also need to enable the `FAST_WAVES` option to tell the compiler that you
specifically want it to include the SIMD version of PlaneWaveTurbulence.
3. Finally, the CPU that will actually run the code needs to support the
required extensions: AVX and everything below. These extensions are relatively
common, but there may still be processors in use that do not support them.

[gj99]: https://doi.org/10.1086/307452
[td13]: https://doi.org/10.1063/1.4789861
 */
class PlaneWaveTurbulence : public TurbulentField {
private:
  int Nm;

  std::vector<Vector3d> xi;
  std::vector<Vector3d> kappa;
  std::vector<double> phi;
  std::vector<double> costheta;
  std::vector<double> beta;
  std::vector<double> Ak;
  std::vector<double> k;

  // data for FAST_WAVES
  int avx_Nm;
  int align_offset;
  std::vector<double> avx_data;
  // the following are index bases into the avx_data array.
  // since each subarray has avx_Nm elements, the start offset
  // of each subarray can be computed by multiplying the two,
  // and then adding on the alignment offset.
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
      Create a new instance of PlaneWaveTurbulence with the specified parameters. This
     generates all of the wavemodes according to the given parameters.
      @param Brms         root mean square field strength for generated field
      @param q, s         the spectral indexes, as defined in the TD13 paper.
     Usually, you'd want to use s=5/3 and q=0 here, which will be comparable to
     passing -11/3 to `initTurbulence`.
      @param lMin, lMax   minimum and maximum wave length
      @param Nm           number of wavemodes that will be used when computing
     the field. A higher value will give a more accurate representation of the
     turbulence, but increase the runtime for getField.
      @param              seed can be used to seed the random number generator
     used to generate the field. This works just like in initTurbulence: a seed
     of 0 will lead to a randomly initialized RNG.
  */
  PlaneWaveTurbulence(const TurbulenceSpectrum &spectrum,
	int Nm = 64, int seed = 0);

  /**
     Evaluates the field at the given position.

     Theoretical runtime is O(Nm), where Nm is the number of wavemodes.
  */
  Vector3d getField(const Vector3d &pos) const;
};

/** @} */

} // namespace crpropa

#endif // CRPROPA_PLANEWAVETURBULENCE_H

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
by [(Giacalone and Jokipii, 1999)][GJ99] and [(Tautz and Dosch, 2013)][TD13].
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
Furthermore, depending on the exact setting, the implementation is somewhat slower
(see Schlegel et al. for details).

 ## Using the SIMD optimization
 In order to mitigate some of the performance impact that is inherent in this
method of field generation, an optimized version is provided. According to our
tests (see the paper above), this version runs 20-30x faster than the baseline
implementation and matches the speed of trilinear interpolation on a grid at
a bit less than 100 wavemodes. To do this, it uses special CPU instructions
called AVX, which are unfortunately not supported by every CPU. On Linux, you
can check whether your CPU supports AVX by using the `lscpu` command, and
searching the flags section for the string "avx". Alternatively, if you are
building on the same CPU that CRPropa will run on, you can have the compiler
detect this for you automatically (see below).

 An additional speedup of about 1.33 can be achieved if the CPU supports the
FMA extension in addition to AVX. Again, you can either check for this
manually, or have the compiler figure it out for you.

 **Note** that the optimized and non-optimized implementations to not return
the exact same results. In fact, since the effective wave numbers used
by the optimized implementation are very slightly different from those
used by the non-optimized version (a difference smaller than the precision
of a double, but nevertheless relevant at some point), the wavemodes go
out of phase for large distances from the origin, and the fields are no longer
comparable at all.

 ### If you are building on the same machine that the code will run on:

1. In cmake: enable the FAST_WAVES flag.
2. Also in cmake: set SIMD_EXTENSIONS to "native". The compiler will automatically
detect support for your CPU and run the build with the appropriate settings.
3. Generate files and exit cmake, then build.
3. If your CPU does not support the necessary extensions, the build will fail
with an error telling you so. In this case, you won't be able to use the optimization;
go back into cmake, disable FAST_WAVES, and build again.
4. If the build runs through without errors, the code is built with the optimization.

 ### If you are building on a different machine from the one where the code will run:

1. Figure out which extensions your target machine supports, using `lscpu | grep avx`
for AVX, and `lscpu | grep fma` for FMA. Run these commands on your target machine,
not on your build machine. If your target machine does not run Linux, you may have to
use different commands.
2. If the CPU does not support AVX, you can stop here, and build normally. The CPU
does not support the necessary extensions, and will not run code that is compiled
using them.
3. If the CPU does support AVX, continue. In cmake, enable the FAST_WAVES flag.
4. Also in cmake, set SIMD_EXTENSIONS to the appropriate value. If your target CPU
supports only AVX, but not FMA, set it to "avx". If it supports both, set it to
"avx+fma".
5. Generate and exit cmake, then build.
6. If you made a mistake and used the wrong flags, you should see an illegal
instruction error when trying to import CRPropa, or (at the latest) when trying
to call `getField`.

[GJ99]: https://doi.org/10.1086/307452
[TD13]: https://doi.org/10.1063/1.4789861
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
	// iAxi is a combined array containing the product of Ak * xi
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
	    Create a new instance of PlaneWaveTurbulence with the specified
	   parameters. This generates all of the wavemodes according to the given
	   parameters.
	    @param spectrum TurbulenceSpectrum
	    @param Nm	number of wavemodes that will be used when computing
	   the field. A higher value will give a more accurate representation of the
	   turbulence, but increase the runtime for getField.
	    @param seed can be used to seed the random number generator
	   used to generate the field. This works just like in initTurbulence: a
	   seed of 0 will lead to a randomly initialized RNG.
	*/
	PlaneWaveTurbulence(const TurbulenceSpectrum &spectrum, int Nm = 64,
	                    int seed = 0);

	/**
	   Evaluates the field at the given position.

	   Theoretical runtime is O(Nm), where Nm is the number of wavemodes.
	*/
	Vector3d getField(const Vector3d &pos) const;
};

/** @} */

} // namespace crpropa

#endif // CRPROPA_PLANEWAVETURBULENCE_H

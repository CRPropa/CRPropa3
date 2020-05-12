#ifndef CRPROPA_SIMPLEGRIDTURBULENCE_H
#define CRPROPA_SIMPLEGRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/Grid.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"

#include "kiss/logger.h"
#include "kiss/string.h"

/**
 @file
 @brief Generate a turbulent field on a grid and related functions.

 This file contains a number of functions related to scalar and vector grids
 (Grid.h).
 */

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class SimpleGridTurbulence
 @brief Turbulent grid-based magnetic field with a simple power-law spectrum
 */
class SimpleGridTurbulence : public TurbulentField {
private:
  double lMin, lMax;
  double boxSize, spacing;
  int gridSize;
  unsigned int seed;
  ref_ptr<Grid3f> gridPtr;

  void initGrid();

public:
  /**
   Create a random initialization of a turbulent field.
   @param Brms	 RMS field strength
   @param sindex	 Spectral index of the energy spectrum in the inertial
   range
   @param lMin	 Minimum physical scale of the turbulence
   @param lMax	 Maximum physical scale of the turbulence
   @param gridSize Grid size (prefer 2^N, where N = int)
   @param boxSize  The length of one side of the box in physical units
   @param seed	 Random seed
   */
  SimpleGridTurbulence(double Brms, double sindex, double lMin, double lMax,
                       int gridSize, double boxSize, unsigned int seed = 0);
  SimpleGridTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin,
                       double lMax, double alpha = -11. / 3.,
                       unsigned int seed = 0);

  /**
      @brief       compute the magnetic field coherence length according to the
     formula in  Harari et Al JHEP03(2002)045
      @return Lc   coherence length of the magnetic field
  */
  double getCorrelationLength() const;
  static double turbulentCorrelationLength(double lMin, double lMax,
                                           double sindex);

  Vector3d getField(const Vector3d &pos) const;

  static void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin,
                             double lMax, double alpha, int seed);
};

// Compatibility with old functions from GridTurbulence:

/** Analytically calculate the correlation length of the simple model turbulent
 * field */
inline double turbulentCorrelationLength(double lMin, double lMax,
                                         double alpha = -11 / 3.) {
  KISS_LOG_WARNING << "turbulentCorrelationLength is deprecated and will be "
                      "removed in the future. Replace it with an appropriate "
                      "turbulent field model and call getCorrelationLength().";
  return SimpleGridTurbulence::turbulentCorrelationLength(lMin, lMax,
                                                          -alpha - 2);
}

/**
 Create a random initialization of a turbulent field.
 @param lMin	Minimum wavelength of the turbulence
 @param lMax	Maximum wavelength of the turbulence
 @param alpha	Power law index of <B^2(k)> ~ k^alpha (alpha = -11/3 corresponds
 to a Kolmogorov spectrum)
 @param Brms	RMS field strength
 @param seed	Random seed
 */
inline void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin,
                           double lMax, double alpha = -11 / 3., int seed = 0) {
  KISS_LOG_WARNING
      << "initTurbulence is deprecated and will be removed in the future. "
         "Replace it with an appropriate turbulent field model instance.";
  SimpleGridTurbulence::initTurbulence(grid, Brms, lMin, lMax, alpha, seed);
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_SIMPLEGRIDTURBULENCE_H

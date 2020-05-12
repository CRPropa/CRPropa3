#ifndef CRPROPA_GRIDTURBULENCE_H
#define CRPROPA_GRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/Grid.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"

#include "fftw3.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class GridTurbulence
 @brief Turbulent grid-based magnetic field with a general energy spectrum
 */
class GridTurbulence : public TurbulentField {
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
   @param qindex	 Spectral index of the energy spectrum in the energy
   range
   @param lBendover Bend-over scale
   @param lMin	 Minimum physical scale of the turbulence
   @param lMax	 Maximum physical scale of the turbulence
   @param gridSize Grid size (prefer 2^N, where N = int)
   @param boxSize  The length of one side of the box in physical units
   @param seed	 Random seed
   */
  GridTurbulence(double Brms, double sindex, double qindex, double lBendover,
                 double lMin, double lMax, int gridSize, double boxSize,
                 unsigned int seed = 0);

  Vector3d getField(const Vector3d &pos) const;

  void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin,
                      double lMax, double alpha, int seed);
};

/* Helper functions for synthetic turbulent field models */

// Check the grid properties before the FFT procedure
void checkGridRequirementsTEMP(ref_ptr<Grid3f> grid, double lMin, double lMax);

// Execute inverse discrete FFT in-place for a 3D grid, from complex to real
// space
void executeInverseFFTInplaceTEMP(ref_ptr<Grid3f> grid, fftwf_complex *Bkx,
                                  fftwf_complex *Bky, fftwf_complex *Bkz);

/**
 Calculate the omnidirectional power spectrum E(k) for a given turbulent field
 Returns a vector of pairs (k_i, E(k_i))
*/
std::vector<std::pair<int, float>> gridPowerSpectrum(ref_ptr<Grid3f> grid);

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_GRIDTURBULENCE_H

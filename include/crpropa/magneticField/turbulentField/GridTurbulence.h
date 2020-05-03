#ifndef CRPROPA_GRIDTURBULENCE_H
#define CRPROPA_GRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/Grid.h"

#include "fftw3.h"

/**
 @file
 @brief Generate a turbulent field on a grid and related functions.

 This file contains a number of functions related to scalar and vector grids (Grid.h).
 */

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class GridTurbulence
 @brief Turbulent grid-based magnetic field with a simple power-law spectrum
 */
class GridTurbulence: public TurbulentField {
private:
	double lMin, lMax;
	double boxSize, spacing;
	int	gridSize;
	unsigned int seed;
	ref_ptr<Grid3f> gridPtr;

	void initGrid();
public:
	GridTurbulence(double Brms,
		double lMin, double lMax, double sindex = 5./3., unsigned int seed = 0);
	GridTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax,
    	double alpha = -11./3., unsigned int seed = 0);

	double getCorrelationLength() const;
	static double turbulentCorrelationLength(double lMin, double lMax, double sindex);

	Vector3d getField(const Vector3d& pos) const;
	
	void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax, double alpha, int seed);
};

/* Helper functions for synthetic turbulent field models */

// Check the grid properties before the FFT procedure
void checkGridRequirementsTEMP(ref_ptr<Grid3f> grid, double lMin, double lMax);

// Execute inverse discrete FFT in-place for a 3D grid, from complex to real space
void executeInverseFFTInplaceTEMP(ref_ptr<Grid3f> grid, fftwf_complex* Bkx, fftwf_complex* Bky, fftwf_complex* Bkz);

/**
 Calculate the omnidirectional power spectrum E(k) for a given turbulent field
 Returns a vector of pairs (k_i, E(k_i))
*/
std::vector<std::pair<int, float> > gridPowerSpectrum(ref_ptr<Grid3f> grid); 

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_GRIDTURBULENCE_H

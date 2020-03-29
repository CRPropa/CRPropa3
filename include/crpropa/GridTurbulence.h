#ifndef CRPROPA_GRIDTURBULENCE_H
#define CRPROPA_GRIDTURBULENCE_H

#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"

/**
 @file
 @brief Generate a turbulent field on a grid and related functions.

 This file contains a number of functions related to scalar and vector grids (Grid.h).
 */

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/** Analytically calculate the correlation length of a turbulent field */
double turbulentCorrelationLength(double lMin, double lMax,
		double alpha = (-11./3.));

#ifdef CRPROPA_HAVE_FFTW3F
/**
 Calculate the omnidirectional power spectrum E(k) for a given turbulent field
 Returns a vector of pairs (k_i, E(k_i))
*/
std::vector<std::pair<int, float> > gridPowerSpectrum(ref_ptr<VectorGrid> grid); 

/**
 Create a random initialization of a turbulent field.
 @param lMin	Minimum wavelength of the turbulence
 @param lMax	Maximum wavelength of the turbulence
 @param alpha	Power law index of <B^2(k)> ~ k^alpha (alpha = -11/3 corresponds to a Kolmogorov spectrum)
 @param Brms	RMS field strength
 @param seed	Random seed
 */
void initTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, 
	   double alpha = -11./3., int seed = 0);

/**
 Same as the normal turbulent field but with helicity.
 @param H	Helicity
*/
void initHelicalTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, 
	   double alpha = -11./3., int seed = 0, double H = 0);

/**
 Same as the normal turbulent field but with the bendover.
 @param lambda	Bendover scale, usually defined as the grid total size divided by number (example: lambda = grid_size/10.)
*/
void initTurbulenceWithBendover(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax,
	   double alpha = -11./3., int seed = 0, double lambda = 1);

#endif // CRPROPA_HAVE_FFTW3F

/** @}*/
} // namespace crpropa

#endif // CRPROPA_GRIDTOOLS_H

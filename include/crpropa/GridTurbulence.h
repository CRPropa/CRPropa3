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
 Create a random initialization of a turbulent field.
 @param lMin	Minimum wavelength of the turbulence
 @param lMax	Maximum wavelength of the turbulence
 @param alpha	Power law index of <B^2(k)> ~ k^alpha (alpha = -11/3 corresponds to a Kolmogorov spectrum)
 @param Brms	RMS field strength
 @param seed	Random seed
 */
void initTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, 
	   double alpha = -11./3., int seed = 0, bool helicity = false, double H = 0, bool bendover = false, double bendover_lambda = 0);

#endif // CRPROPA_HAVE_FFTW3F

/** @}*/
} // namespace crpropa

#endif // CRPROPA_GRIDTOOLS_H

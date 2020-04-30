#ifndef CRPROPA_GRIDTURBULENCE_H
#define CRPROPA_GRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"

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
 @brief Provides a turbulent magnetic field generated on a grid
 */
class GridTurbulence: public TurbulentField {
private:
	double lMin, lMax;
	MagneticFieldGrid mfield;

	initTurbulence();
public:

	Vector3d getField(const Vector3d& pos) const;
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_GRIDTOOLS_H

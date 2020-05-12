#ifndef CRPROPA_SIMPLEGRIDTURBULENCE_H
#define CRPROPA_SIMPLEGRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/Grid.h"

#include "kiss/string.h"
#include "kiss/logger.h"

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
 @class SimpleGridTurbulence
 @brief Turbulent grid-based magnetic field with a simple power-law spectrum
 */
class SimpleGridTurbulence: public TurbulentField {
private:
	double lMin, lMax;
	double boxSize, spacing;
	int	gridSize;
	unsigned int seed;
	ref_ptr<Grid3f> gridPtr;

	void initGrid();
public:
	SimpleGridTurbulence(double Brms,
		double lMin, double lMax, double sindex = 5./3., unsigned int seed = 0);
	SimpleGridTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax,
    	double alpha = -11./3., unsigned int seed = 0);

	double getCorrelationLength() const;
	static double turbulentCorrelationLength(double lMin, double lMax, double sindex);

	Vector3d getField(const Vector3d& pos) const;
	
	static void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax, double alpha, int seed);
};

// Compatibility with old functions from GridTurbulence:

/** Analytically calculate the correlation length of the simple model turbulent field */
inline double turbulentCorrelationLength(double lMin, double lMax, double alpha = -11/3.) {
	KISS_LOG_WARNING << "turbulentCorrelationLength is deprecated and will be removed in the future. Replace it with an appropriate turbulent field model and call getCorrelationLength().";
	return SimpleGridTurbulence::turbulentCorrelationLength(lMin, lMax, -alpha - 2);
}

inline void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax,
		double alpha = -11/3., int seed = 0) {
	KISS_LOG_WARNING << "initTurbulence is deprecated and will be removed in the future. Replace it with an appropriate turbulent field model instance.";
	SimpleGridTurbulence::initTurbulence(grid, Brms, lMin, lMax, alpha, seed);
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_SIMPLEGRIDTURBULENCE_H

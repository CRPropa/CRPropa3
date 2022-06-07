#ifndef CRPROPA_SIMPLEGRIDTURBULENCE_H
#define CRPROPA_SIMPLEGRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/magneticField/turbulentField/GridTurbulence.h"

#include "kiss/logger.h"
#include "kiss/string.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class SimpleTurbulenceSpectrum
 @brief Defines the energy spectrum of simple power-law turbulence
 */
class SimpleTurbulenceSpectrum : public TurbulenceSpectrum {
	const int constScaleBendover = 1000; 	// define the bandover scale as 1000 * lMax to ensure k * lBendover >> 1. The bendover scale is necessary for the implementation of PlaneWaveTurbulence. 
  public:
	/**
	 @param Brms		root mean square field strength for generated field
	 @param lMin	 	Minimum physical scale of the turbulence
	 @param lMax	 	Maximum physical scale of the turbulence
	 @param lBendover	the bend-over scale is set to 1000 times lMax to ensure to be in the inertial range. This should not be changed.
	 @param sindex	 	Spectral index of the energy spectrum in the inertial range
	*/
	SimpleTurbulenceSpectrum(double Brms, double lMin, double lMax,
	                         double sIndex = 5. / 3)
	    : TurbulenceSpectrum(Brms, lMin, lMax, constScaleBendover * lMax, sIndex, 0) {}
	~SimpleTurbulenceSpectrum() {}

	/**
	General energy spectrum for synthetic turbulence models
	*/
	double energySpectrum(double k) const {
		return std::pow(k, -getSindex() - 2);
	}

	/**
	    @brief       compute the magnetic field coherence length according to
	   the formula in  Harari et al. JHEP03(2002)045
	    @return Lc   coherence length of the magnetic field
	*/
	double getCorrelationLength() const {
		return turbulentCorrelationLength(getLmin(), getLmax(),
		                                  getSindex());
	}
	static double turbulentCorrelationLength(double lMin, double lMax,
	                                         double s) {
		double r = lMin / lMax;
		return lMax / 2 * (s - 1) / s * (1 - pow(r, s)) / (1 - pow(r, s - 1));
	}
};

/**
 @class SimpleGridTurbulence
 @brief Turbulent grid-based magnetic field with a simple power-law spectrum
 */
class SimpleGridTurbulence : public GridTurbulence {
  public:
	/**
	 Create a random initialization of a turbulent field.
	 @param spectrum    TurbulenceSpectrum instance to define the spectrum of
	 turbulence
	 @param gridProp	GridProperties instance to define the underlying grid
	 @param seed	 Random seed
	 */
	SimpleGridTurbulence(const SimpleTurbulenceSpectrum &spectrum,
	                     const GridProperties &gridProp, unsigned int seed = 0);

	static void initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin,
	                           double lMax, double alpha, int seed);
};

// Compatibility with old functions from GridTurbulence:

/** Analytically calculate the correlation length of the simple model turbulent
 * field */
inline double turbulentCorrelationLength(double lMin, double lMax,
                                         double alpha = -11 / 3.) {
	KISS_LOG_WARNING
	    << "turbulentCorrelationLength is deprecated and will be "
	       "removed in the future. Replace it with a more appropriate "
	       "turbulent field model and call getCorrelationLength().";
	return SimpleTurbulenceSpectrum::turbulentCorrelationLength(lMin, lMax,
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
	       "Replace it with a more appropriate turbulent field model instance.";
	SimpleGridTurbulence::initTurbulence(grid, Brms, lMin, lMax, alpha, seed);
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_SIMPLEGRIDTURBULENCE_H

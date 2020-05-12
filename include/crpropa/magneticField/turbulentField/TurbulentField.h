#ifndef CRPROPA_TURBULENTFIELD_H
#define CRPROPA_TURBULENTFIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include <cmath>

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class TurbulenceSpectrum
 @brief Defines the energy spectrum of turbulence
 */
class TurbulenceSpectrum : public Referenced {
  private:
	const double Brms; /**< Brms value of the turbulent field (normalization) */
	const double sIndex; /**< Spectral index for the inertial range, for example
	                  s=5/3 for Kolmogorov spectrum */
	const double qIndex; /**< Spectral index for the injection range, for
	                  example q=4 for 3D homogeneous turbulence */
	const double lBendover;  /**< the bend-over scale */
	const double lMin, lMax; /**< Min and Max scale of turbulence */

  public:
	/**
	   @param Brms         root mean square field strength for generated field
	   @param lMin, lMax   minimum and maximum wave length
	   @param lBendover	   the bend-over scale
	   @param q, s         the spectral indices
	                       Usually, you'd want to use s=5./3.
	                       and q=4 here, which will be comparable
	*/
	TurbulenceSpectrum(double Brms, double lMin, double lMax,
	                   double lBendover = 1, double sIndex = 5. / 3.,
	                   double qIndex = 4)
	    : Brms(Brms), lMin(lMin), lMax(lMax), lBendover(lBendover),
	      sIndex(sIndex), qIndex(qIndex) {
		if (lMin > lMax) {
			throw std::runtime_error("TurbulenceSpectrum: lMin > lMax");
		}
		if (lMin <= 0) {
			throw std::runtime_error("TurbulenceSpectrum: lMin <= 0");
		}
	}

	~TurbulenceSpectrum() {}

	double getBrms() const { return Brms; }
	double getLmin() const { return lMin; }
	double getLmax() const { return lMax; }
	double getLbendover() const { return lBendover; }
	double getSindex() const { return sIndex; }
	double getQindex() const { return qIndex; }

	/**
	General energy spectrum for synthetic turbulence models
	*/
	virtual double energySpectrum(double k) const {
		return lBendover * std::pow(k * lBendover, qIndex) /
		       std::pow(1.0 + k * k * lBendover * lBendover,
		                (sIndex + qIndex) / 2);
	}

	/**
	Normalization for the above defined energy spectrum
	*/
	double spectrumNormalization() const {
		return std::tgamma((sIndex + qIndex) / 2.0) /
		       (2.0 * std::tgamma((sIndex - 1) / 2.0) *
		        std::tgamma((qIndex + 1) / 2.0));
	}

	/**
  Computes the magnetic field coherence length
	Obtained from the definition of \f$l_c = 1/B_{\rm rms}^2 \int_0^\infty dr
  \langleB(0)B^*(r)\rangle$
	Approximates the true value correctly as long as lBendover <= lMax/8 (~5%
  error) (for the true value the above integral should go from lMin to lMax)
	*/
	virtual double getCorrelationLength() const {
		return 4 * M_PI / ((sIndex + 2.0) * sIndex) * spectrumNormalization() *
		       lBendover;
	}
};

/**
 @class TurbulentField
 @brief An abstract base class for different models of turbulent magnetic fields

 This module provides common methods for all turbulent (synthetic) magnetic
 fields. Does not actually implement any turbulent field.
 */
class TurbulentField : public MagneticField {
  protected:
	const TurbulenceSpectrum &spectrum;

  public:
	TurbulentField(const TurbulenceSpectrum &spectrum) : spectrum(spectrum) {}
	virtual ~TurbulentField() {}
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_TURBULENTFIELD_H

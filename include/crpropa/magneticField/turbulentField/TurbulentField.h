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
 @brief Defines the energy spectrum of turbulence parametrizied by A(k) ~ k^q /(1 + k^2)^{(s + q)/2 + 1}
 */
class TurbulenceSpectrum : public Referenced {
  private:
	const double Brms; /**< Brms value of the turbulent field (normalization) */
	const double sIndex; /**< Spectral index for the inertial range, for example
	                  s=5/3 for Kolmogorov spectrum; in some parts of the code this
			  parameter is referred by alpha which is the total 3D isotropic
			  spectrum with additional k^2 and the minus sign, e.g.,
			  for Kolmogorov: alpha = -(s + 2) */
	const double qIndex; /**< Spectral index for the injection range, for
	                  example q=4 for 3D homogeneous turbulence */
	const double lBendover;  /**< the bend-over scale */
	const double lMin, lMax; /**< Min and Max scale of turbulence */

  protected:
	/**
	Normalization for the below defined Lc
	*/
	double spectrumNormalization() const {
		return std::tgamma((sIndex + qIndex) / 2.0) /
		       (2.0 * std::tgamma((sIndex - 1) / 2.0) *
		        std::tgamma((qIndex + 1) / 2.0));
	}

  public:
	/**
	   @param Brms         root mean square field strength for generated field
	 @param lMin	 Minimum physical scale of the turbulence
	 @param lMax	 Maximum physical scale of the turbulence
	   @param lBendover	   the bend-over scale
	 @param sindex	 Spectral index of the energy spectrum in the inertial
	 range
	 @param qindex	 Spectral index of the energy spectrum in the energy
	 range
	*/
	TurbulenceSpectrum(double Brms, double lMin, double lMax,
	                   double lBendover = 1, double sIndex = (5. / 3.),
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
	General energy spectrum for synthetic turbulence models (not normalized!)
	with normalized ^k = k*lBendover
	*/
	virtual double energySpectrum(double k) const {
		double kHat = k * lBendover;
		return std::pow(kHat, qIndex) /
				       std::pow(1.0 + kHat * kHat,
			                (sIndex + qIndex) / 2.0 + 1.0);
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

	double getBrms() const { return spectrum.getBrms(); }
	virtual double getCorrelationLength() const {
		return spectrum.getCorrelationLength();
	}
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_TURBULENTFIELD_H

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
 @class TurbulentField
 @brief An abstract base class for different models of turbulent magnetic fields

 This module provides common methods for all turbulent (synthetic) magnetic
 fields. Does not actually implement any turbulent field.
 */
class TurbulentField : public MagneticField {
protected:
  double Brms;   /**< Brms value of the turbulent field (normalization) */
  double sindex; /**< Spectral index for the inertial range, for example s=5/3
                    for Kolmogorov spectrum */
  double qindex; /**< Spectral index for the injection range, for example q=4
                    for 3D homogeneous turbulence */
  double lBendover; /**< the bend-over scale */

  /**
  General energy spectrum for synthetic turbulence models
  */
  virtual double energySpectrum(double k) const {
    return lBendover * std::pow(k * lBendover, qindex) /
           std::pow(1.0 + k * k * lBendover * lBendover,
                    (sindex + qindex) / 2);
  }

  /**
  Normalization for the above defined energy spectrum
  */
  double spectrumNormalization() const {
    return std::tgamma((sindex + qindex) / 2.0) /
           (2.0 * std::tgamma((sindex - 1) / 2.0) *
            std::tgamma((qindex + 1) / 2.0));
  }

public:
  TurbulentField(double Brms, double sindex, double qindex = 4,
                 double lBendover = 1)
      : Brms(Brms), sindex(sindex), qindex(qindex), lBendover(lBendover) {
  }

  virtual ~TurbulentField() {}

  double getBrms() const { return Brms; }

  /**
Computes the magnetic field coherence length
  Obtained from the definition of \f$l_c = 1/B_{\rm rms}^2 \int_0^\infty dr
\langleB(0)B^*(r)\rangle$
  Approximates the true value correctly as long as lBendover <= lMax/8 (~5% error)
  (for the true value the above integral should go from lMin to lMax)
  */
  virtual double getCorrelationLength() const {
    return 4 * M_PI / ((sindex + 2.0) * sindex) * spectrumNormalization() *
           lBendover;
  }
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_TURBULENTFIELD_H

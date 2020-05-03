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

 This module provides common methods for all turbulent (synthetic) magnetic fields
 */
class TurbulentField: public MagneticField {
protected:
	double Brms;
	double sindex, qindex;
	double l_bendover;

	virtual double energySpectrum(double k) const {
		return l_bendover * std::pow(k * l_bendover, qindex) /
					std::pow(1.0 + k*k*l_bendover*l_bendover, (sindex+qindex)/2);
	}
	double spectrumNormalization() const {
		return std::tgamma((sindex+qindex)/2.0) /
			(2.0 * std::tgamma((sindex-1)/2.0) * std::tgamma((qindex+1)/2.0));
	}
public:
	TurbulentField(double Brms_, double sindex_, double qindex_ = 4, double l_bendover_ = 1)
		: Brms(Brms_), sindex(sindex_), qindex(qindex_), l_bendover(l_bendover_) {
	}

    virtual ~TurbulentField() {
    }

	double getBrms() const {
		return Brms;
	}

	virtual double getCorrelationLength() const {
		return 4*M_PI/((sindex+2.0)*sindex) * spectrumNormalization() * l_bendover;
	}
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_TURBULENTFIELD_H

#ifndef CRPROPA_TURBULENTFIELD_H
#define CRPROPA_TURBULENTFIELD_H

#include "crpropa/magneticField/MagneticField.h"

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
	virtual double energySpectrum(double s,
								  double q,
								  double l_bendover,
								  double k) const {
		return l_bendover * std::pow(k * l_bendover, q) /
					std::pow(1.0 + k*k*l_bendover*l_bendover, (s+q)/2);
	}
public:
	TurbulentField() {
	}
    virtual ~TurbulentField() {
    }
	
	virtual double getCorrelationLength() const = 0;
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_TURBULENTFIELD_H

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
 @brief An abstract class for different models of turbulent magnetic fields

 This module provides common methods for all turbulent (synthetic) magnetic fields
 */
class TurbulentField: public MagneticField {
public:
	TurbulentField() {};
	virtual double getCorrelationLength() const = 0;
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_TURBULENTFIELD_H

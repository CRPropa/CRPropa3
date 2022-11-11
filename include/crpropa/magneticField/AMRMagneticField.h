#ifndef CRPROPA_AMRMAGNETICFIELD_H
#define CRPROPA_AMRMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_SAGA

#include <iostream>
#include <string>
#include <cstdio>

#ifdef _OPENMP
	#include "omp.h"
#endif

#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Vector3.h"

#include "kiss/logger.h"

#include "saga/LocalProperties.h"
#include "saga/AMRgrid.h"
#include "saga/MagneticField.h"
#include "saga/Referenced.h"



namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class AMRMagneticField
 @brief Wrapper for saga::MagneticField
 Can be found in github.com/rafaelab/saga

 Deprecation Warning: 
 As SAGA (SQLite AMR Grid Application) is no longer supported, the AMRMagenticField
 class will be removed in the future.
 */
class AMRMagneticField: public MagneticField {

private:
	saga::ref_ptr<saga::MagneticField> field;
	double cfLength;
	double cfDensity;
	double cfMagneticField;

public:
	/** Constructor
	@param field_	saga magnetic field 
	@param convLength	rescaling length 
	@param convDensity	rescaling density
	@param convMagneticField rescaling magnetic field strength
*/
	AMRMagneticField(saga::ref_ptr<saga::MagneticField> field_, double convLength, double convDensity, double convMagneticField)
	{
		field = field_;
		cfLength = convLength;
		cfDensity = convDensity;
		cfMagneticField = convMagneticField;

		KISS_LOG_WARNING << "DEPRECATION WARNING: AMRMagneticField class will be removed in the future,"
		<<"as the underlying library (saga) is no longer supported.";
	}

	Vector3d getField(const Vector3d &position) const {

		double x = position.x / cfLength;
		double y = position.y / cfLength;
		double z = position.z / cfLength;

		std::vector<double> b = field->getField(x, y, z);
		Vector3d B;
		B.setXYZ(b[0], b[1], b[2]);
		B = B * cfMagneticField;

		return B;
	}

};
/** @} */

} // namespace crpropa

#endif // CRPROPA_HAVE_SAGA
#endif // CRPROPA_AMRMAGNETICFIELD_H

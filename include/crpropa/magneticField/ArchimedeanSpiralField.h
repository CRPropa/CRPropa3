#ifndef CRPROPA_ARCHIMEDEANSPIRALFIELD_H
#define CRPROPA_ARCHIMEDEANSPIRALFIELD_H

#include "crpropa/magneticField/MagneticField.h"


#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"

namespace crpropa {

/**
 * \addtogroup MagneticFields
 * @{
 */

/**

@class ArchimedeanSpiralField
@brief Magnetic field model following an Archimedean spiral.

See e.g. Jokipii, Levy & Hubbard 1977
*/

class ArchimedeanSpiralField: public MagneticField {
private:
	double B_0; // Magnetic field strength in radial direction at R_0
	double R_0; // Reference level
	double Omega; // Angular velocity of the rotation
	double V_w; // Asymptotic wind speed

public:
/** Constructor
	@param B_0	Magnetic field strength in radial direction at R_0
	@param R_0	Reference level
	@param Omega	Angular velocity of the rotation
	@param V_w	Asymptotic wind speed
*/
	ArchimedeanSpiralField(double B_0, double R_0, double Omega, double V_w);

	Vector3d getField(const Vector3d &pos) const;	
		
	void setB0(double B);
	void setR0(double R);
	void setOmega(double Om);
	void setVw(double v);

	double getB0() const;
	double getR0() const;
	double getOmega() const;
	double getVw() const;
};
/** @} */
	 
} // end namespace crpropa

#endif // CRPROPA_ACHIMEDEANSPIRALFIELD_H

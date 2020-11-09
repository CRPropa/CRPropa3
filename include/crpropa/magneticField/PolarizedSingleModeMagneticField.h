#ifndef CRPROPA_POLARIZEDSINGLEMODEMAGNETICFIELD_H
#define CRPROPA_POLARIZEDSINGLEMODEMAGNETICFIELD_H

#include "crpropa/magneticField/MagneticField.h"

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"

#include <stdexcept>

namespace crpropa {

/**

@class PolarizedSingleModeMagneticField
@brief General poralized single mode magnetic field (including linear, circular and elliptic polarizations)


*/

class PolarizedSingleModeMagneticField: public MagneticField {
private:
	double B_0; // Magnetic field strength in the direction of e_1 at r_0
	double wavelength; // correlation length
	Vector3d r_0; // Reference position
	Vector3d e_1; // First vector spanning the polarization plane
	Vector3d e_2; // Second vector spanning the polarization plane 

public:
/** Constructor
        @param B_0 // Magnetic field strength in the direction of e_1 at r_0
        @param wavelength; // correlation length
        @param r_0; // reference position
        @param e_1; // First vector spanning the polarization plane 
        @param e_2; // Second vector spanning the polarization plane
*/
	PolarizedSingleModeMagneticField(double B_0, double wavelength, Vector3d r_0, Vector3d e_1, Vector3d e_2);

	Vector3d getField(const Vector3d &pos, const double &sigma) const;

	Vector3d getBrmsMaxHelFracField(const Vector3d &pos, const double &fH) const;
	Vector3d getOrthogonalBrmsMaxHelFracField(const Vector3d &pos, const double &fH) const;

	Vector3d getGeneralOrthogonalEllipticSingleModeB0MagneticField(const Vector3d &pos, const double &sigma) const;
	Vector3d getGeneralOrthogonalEllipticSingleModeBrmsMagneticField(const Vector3d &pos, const double &sigma) const;

	Vector3d getSpecialOrthogonalEllipticSingleModeB0MagneticField(const Vector3d &pos, const double &sigma) const;
	Vector3d getSpecialOrthogonalEllipticSingleModeBrmsMagneticField(const Vector3d &pos, const double &sigma) const;

	Vector3d getCircularSingleModeMagneticField(const Vector3d &pos, const double &sigma) const;

	Vector3d getLinearSingleModeB0MagneticField(const Vector3d &pos) const;	
	Vector3d getLinearSingleModeBrmsMagneticField(const Vector3d &pos) const;

	void setB0(double B);
	void setWavelength(double wavelength);
	void setr0(Vector3d r0);
	void sete1(Vector3d e1);
	void sete2(Vector3d e2);

	double getB0() const;
	double getWavelength() const;
	Vector3d getr0() const;
	Vector3d gete1() const;
	Vector3d gete2() const;
};

	 
} // end namespace crpropa

#endif // CRPROPA_POLARIZEDSINGLEMODEMAGNETICFIELD_H

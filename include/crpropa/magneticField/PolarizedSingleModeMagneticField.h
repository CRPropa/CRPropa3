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

namespace crpropa {

/**

@class PolarizedSingleModeMagneticField
@brief 	

*/

/**
	@class PolarizedSingleModeMagneticField
	@brief General poralized single mode magnetic field (including linear, circular and elliptic polarizations and the case of maximal helicity).
 */
  
class PolarizedSingleModeMagneticField: public MagneticField {

private:

	double B_0; // Magnetic field strength in the direction of e_1 at r_0
	double wavelength; // Wavelength of the single mode (corresponds to its coherence length)
	double sigma; // Polarization parameter
	Vector3d r_0; // Reference position
	Vector3d e_1; // First vector spanning the polarization plane
	Vector3d e_2; // Second vector spanning the polarization plane
	std::string flagAmplitudeRms; // Flag to specify whether B_0 denotes the maximum ("amplitude") or the RMS ("rms") value of the magnetic field.
	std::string flagPolarizationHelicity; // Flag to specify whether sigma denotes the standard polarization parameter ("polarization") or f_H, the fraction of maximal helicity.
	std::string flagMode; // Flag to specify the polarization mode; possible choices are "elliptical", "circular" or "linear".

public: 

	PolarizedSingleModeMagneticField(const double &B_0, const double &wavelength, const double &sigma, const Vector3d &r_0, const Vector3d &e_1, const Vector3d &e_2, std::string flagAmplitudeRms, std::string flagPolarizationHelicity, std::string flagMode ) :
			B_0(B_0), wavelength(wavelength), sigma(sigma), r_0(r_0), e_1(e_1), e_2(e_2), flagAmplitudeRms(flagAmplitudeRms), flagPolarizationHelicity(flagPolarizationHelicity), flagMode(flagMode) {
	}

        Vector3d getField(const Vector3d &position) const;
};

} // end namespace crpropa

#endif // CRPROPA_POLARIZEDSINGLEMODEMAGNETICFIELD_H

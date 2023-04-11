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
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class PolarizedSingleModeMagneticField
 @brief General polarized single mode magnetic field (including linear, circular and elliptic polarizations and the case of maximal helicity).
 */
class PolarizedSingleModeMagneticField: public MagneticField {
private:
	//quantities provided by the user:
	double B_0; // Magnetic field strength in the direction of e_1 at r_0 (for flagAmplitudeRms = "amplitude"), or the RMS value of the magnetic field (for flagAmplitudeRms = "rms")
	double wavelength; // Wavelength of the single mode (corresponds to its coherence length)
	double sigma; // Polarization parameter
	Vector3d r_0; // Reference position
	Vector3d e_1; // First vector spanning the polarization plane
	Vector3d e_2; // Second vector spanning the polarization plane
	std::string flagAmplitudeRms; // Flag to specify whether B_0 denotes the maximum ("amplitude") or the RMS ("rms") value of the magnetic field
	std::string flagPolarizationHelicity; // Flag to specify whether sigma denotes the standard polarization parameter ("polarization") or f_H, the fraction of maximal helicity ("helicity")
	std::string flagMode; // Flag to specify the polarization mode; possible choices are "elliptical", "circular" or "linear"

	//derived quantities:
	Vector3d unitVector_1; // Normalized vector e_1
	Vector3d unitVector_2; // Normalized vector e_2
	Vector3d wavevector; // Wavevector of the mode (proportional to e_2 cross e_1)
	double B_max; // Maximal value of the magnetic field (i.e. the amplitude/semi-major value of the mode)

public:
	/**
	 * Constructor
	 * @param B_0 						Magnetic field strength in the direction of e_1 at r_0 (for flagAmplitudeRms = "amplitude"), or the RMS value of the magnetic field (for flagAmplitudeRms = "rms")	
	 * @param wavelength				Wavelength of the single mode (corresponds to its coherence length)
	 * @param sigma 					Polarization parameter
	 * @param r_0						Reference position
	 * @param e_1						First vector spanning the polarization plane
	 * @param e_2	 					Second vector spanning the polarization plane
	 * @param flagAmplitudeRms			Flag to specify whether B_0 denotes the maximum ("amplitude") or the RMS ("rms") value of the magnetic field
	 * @param flagPolarizationHelicity	Flag to specify whether sigma denotes the standard polarization parameter ("polarization") or f_H, the fraction of maximal helicity ("helicity")
	 * @param flagMode 					Flag to specify the polarization mode; possible choices are "elliptical", "circular" or "linear"
	*/
	PolarizedSingleModeMagneticField(const double &B_0, const double &wavelength, const double &sigma, const Vector3d &r_0, const Vector3d &e_1, const Vector3d &e_2, std::string flagAmplitudeRms, std::string flagPolarizationHelicity, std::string flagMode);

	Vector3d getField(const Vector3d &position) const;
};
/** @} */

} // end namespace crpropa

#endif // CRPROPA_POLARIZEDSINGLEMODEMAGNETICFIELD_H

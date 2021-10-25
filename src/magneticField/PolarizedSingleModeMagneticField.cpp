#include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"

namespace crpropa {

Vector3d PolarizedSingleModeMagneticField::getField(const Vector3d &position) const {
	if (flagMode == "elliptical") {
		if (abs(sigma) > 1)
			throw std::runtime_error("PolarizedSingleModeMagneticField: The value of the  polarization parameter has to lie in the range [-1;+1].");
	}
	else if (flagMode == "circular") {
		if (abs(sigma) != 1)
			throw std::runtime_error("PolarizedSingleModeMagneticField: For circular polarization the value of the polarization parameter has to be equal to -1 or +1.");
	}
	else if (flagMode == "linear") {
		if (abs(sigma) != 0)
			throw std::runtime_error("PolarizedSingleModeMagneticField: For linear polarization the value of the polarization parameter has to be equal to 0.");
	}
	else
		throw std::runtime_error("PolarizedSingleModeMagneticField: Wrong value for flagMode. Please choose \"elliptical\" or \"circular\" or \"linear\".");

	if (e_1.dot(e_2) != 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: e_1 and e_2 have to be orthogonal to each other.");

	Vector3d unitVector_1;
	if (e_1.getR() == 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: Vector e_1 cannot be zero.");
	unitVector_1 = e_1 / e_1.getR();

	Vector3d unitVector_2;
	if (e_2.getR() == 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: Vector e_2 cannot be zero.");
	unitVector_2 = e_2 / e_2.getR();

	Vector3d wavevector = e_2.cross(e_1);

	//This check is necessary if a polarization with non-orthogonal spanning vectors is desired. This is not implemented in this version, so any corresponding modifications should be made with caution.
	//if (wavevector.getR() == 0) 
	//	throw std::runtime_error("PolarizedSingleModeMagneticField: e_1 cannot be parallel to e_2.");

	wavevector = wavevector / wavevector.getR();

	if (wavelength == 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: The correlation length cannot be zero.");
	wavevector = wavevector/wavevector.getR() * 2 * M_PI / wavelength;

	if (!((flagPolarizationHelicity == "helicity") || (flagPolarizationHelicity == "polarization")))
		throw std::runtime_error("PolarizedSingleModeMagneticField: Wrong value for flagPolarizationHelicity. Please choose \"polarization\" or \"helicity\".");
	if (flagPolarizationHelicity == "helicity" && abs(sigma) != 1)
		throw std::runtime_error("PolarizedSingleModeMagneticField: In helicity mode only the maximum helicity case (sigma = +1 or sigma = -1) may be chosen.");

	double prefactor;

	if (flagAmplitudeRms == "amplitude")
		prefactor = 1;
	else if (flagAmplitudeRms == "rms")
		prefactor = sqrt(2 / (1 + sigma * sigma));
	else
		throw std::runtime_error("PolarizedSingleModeMagneticField: Wrong value for flagAmplitudeRms. Please choose \"amplitude\" or \"rms\".");

	Vector3d delta_r = position - r_0;
  
  return prefactor * B_0 * (unitVector_1 * cos(wavevector.dot(delta_r)) + unitVector_2 * sigma * sin(wavevector.dot(delta_r)));
}

} //end namespace crpropa

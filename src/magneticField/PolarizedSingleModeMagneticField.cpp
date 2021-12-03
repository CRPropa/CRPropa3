#include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"

namespace crpropa {

PolarizedSingleModeMagneticField::PolarizedSingleModeMagneticField( const double &B_0, const double &wavelength, const double &sigma, const Vector3d &r_0, const Vector3d &e_1, const Vector3d &e_2, std::string flagAmplitudeRms, std::string flagPolarizationHelicity, std::string flagMode ) {
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

	if (e_1.getR() == 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: Vector e_1 cannot be zero.");
	unitVector_1 = e_1 / e_1.getR();

	if (e_2.getR() == 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: Vector e_2 cannot be zero.");
	unitVector_2 = e_2 / e_2.getR();

	wavevector = e_2.cross(e_1);

	//This check is necessary if a polarization with non-orthogonal spanning vectors is desired. This is not implemented in this version, so any corresponding modifications should be made with caution.
	//if (wavevector.getR() == 0)
	//        throw std::runtime_error("PolarizedSingleModeMagneticField: e_1 cannot be parallel to e_2.");

	wavevector = wavevector / wavevector.getR();

	if (wavelength == 0)
		throw std::runtime_error("PolarizedSingleModeMagneticField: The correlation length cannot be zero.");
	wavevector = wavevector * 2 * M_PI / wavelength;

	if (!((flagPolarizationHelicity == "helicity") || (flagPolarizationHelicity == "polarization")))
		throw std::runtime_error("PolarizedSingleModeMagneticField: Wrong value for flagPolarizationHelicity. Please choose \"polarization\" or \"helicity\".");
	if (flagPolarizationHelicity == "helicity" && abs(sigma) != 1)
		throw std::runtime_error("PolarizedSingleModeMagneticField: In helicity mode only the maximum helicity case (sigma = +1 or sigma = -1) may be chosen.");

	if (flagAmplitudeRms == "amplitude")
		B_max = B_0;
	else if (flagAmplitudeRms == "rms")
		B_max = sqrt(2 / (1 + sigma * sigma)) * B_0;
	else
		throw std::runtime_error("PolarizedSingleModeMagneticField: Wrong value for flagAmplitudeRms. Please choose \"amplitude\" or \"rms\".");

	this->r_0 = r_0;
	this->sigma = sigma;
}

Vector3d PolarizedSingleModeMagneticField::getField(const Vector3d &position) const {
	Vector3d delta_r = position - r_0;

	return B_max * (unitVector_1 * cos(wavevector.dot(delta_r)) + unitVector_2 * sigma * sin(wavevector.dot(delta_r)));
}

} //end namespace crpropa

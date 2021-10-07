#include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"

namespace crpropa {

PolarizedSingleModeMagneticField::PolarizedSingleModeMagneticField(double B_0, double wavelength, Vector3d r_0, Vector3d e_1, Vector3d e_2) {
	setB0(B_0);
	setWavelength(wavelength);
	setr0(r_0);
	sete1(e_1);
	sete2(e_2);
}

Vector3d PolarizedSingleModeMagneticField::getField(const Vector3d &pos, const double &sigma) const {
	Vector3d deltar = pos - r_0;

	Vector3d unitVector_1;
	if (e_1.getR() == 0)
		throw std::runtime_error("Vector e_1 cannot be zero.");
	unitVector_1 = e_1 / e_1.getR();

	Vector3d unitVector_2;
	if (e_2.getR() == 0)
		throw std::runtime_error("Vector e_2 cannot be zero.");
	unitVector_2 = e_2 / e_2.getR();

	Vector3d wavevector = e_2.cross(e_1);
	if (wavevector.getR() == 0)
		throw std::runtime_error("e_1 cannot be parallel to e_2.");
	wavevector = wavevector / wavevector.getR();

	if (wavelength == 0)
		throw std::runtime_error("The correlation length cannot be zero.");
	wavevector = wavevector/wavevector.getR() * 2 * M_PI / wavelength;

	return B_0 * (unitVector_1 * cos(wavevector.dot(deltar)) + unitVector_2 * sigma * sin(wavevector.dot(deltar)));
}

Vector3d PolarizedSingleModeMagneticField::getMaxHelFracField(const Vector3d &pos, const double &fH, const std::string &Bflag) const {
	if (abs(fH) > 1)
		throw std::runtime_error("The value of sigma has to be in the range [-1;1].");

	if (Bflag == "amplitude") {
		getField(pos, fH);
	}
	else if (Bflag == "rms") {
		if (fH == 0) {
			getField(pos, 0);
		}
		else {
			return sqrt(1 + sqrt(1 - fH * fH))*getField(pos,(1 - sqrt(1 - fH * fH))/fH);
		}
	}
	else {
		throw std::runtime_error("Wrong value for Bflag. Please choose \"B0\" or \"Brms\".");
	}
}

Vector3d PolarizedSingleModeMagneticField::getOrthogonalMaxHelFracField(const Vector3d &pos, const double &fH, const std::string &Bflag) const {
	if (e_1.dot(e_2) != 0)
		throw std::runtime_error("e_1 and e_2 have to be orthogonal to each other.");

	return getMaxHelFracField(pos, fH, Bflag);	
}

Vector3d PolarizedSingleModeMagneticField::getGeneralOrthogonalEllipticField(const Vector3d &pos, const double &sigma, const std::string &Bflag) const {
	if (e_1.dot(e_2) != 0)
		throw std::runtime_error("e_1 and e_2 have to be orthogonal to each other.");

	if (Bflag == "amplitude") {
		return getField(pos, sigma);
	}
	else if (Bflag == "rms") {
		return sqrt(2 / (1 + sigma * sigma)) * getField(pos, sigma);
	}
	else {
		throw std::runtime_error("Wrong value for Bflag. Please choose \"B0\" or \"Brms\".");
	}
}

Vector3d PolarizedSingleModeMagneticField::getSpecialOrthogonalEllipticField(const Vector3d &pos, const double &sigma, const std::string &Bflag) const {
	if (abs(sigma) > 1)
		throw std::runtime_error("The value of sigma has to be in the range [-1;1].");

	return getGeneralOrthogonalEllipticField(pos, sigma, Bflag);
}

Vector3d PolarizedSingleModeMagneticField::getOrthogonalCircularField(const Vector3d &pos, const double &sigma) const {
	if (sigma == 0)
		throw std::runtime_error("The value of sigma cannot be zero.");

	return getGeneralOrthogonalEllipticField(pos, sigma/abs(sigma), "B0");
}

Vector3d PolarizedSingleModeMagneticField::getLinearField(const Vector3d &pos, const std::string &Bflag) const {
	return getGeneralOrthogonalEllipticField(pos, 0, Bflag);
}


void PolarizedSingleModeMagneticField::setB0(double B) {
	B_0 = B;
	return;
}

void PolarizedSingleModeMagneticField::setWavelength(double wavelen) {
	wavelength = wavelen;
	return;
}

void PolarizedSingleModeMagneticField::setr0(Vector3d r0) {
	r_0 = r0;
	return;
}

void PolarizedSingleModeMagneticField::sete1(Vector3d e1) {
	e_1 = e1;
	return;
}

void PolarizedSingleModeMagneticField::sete2(Vector3d e2) {
	e_2 = e2;
	return;
}

}

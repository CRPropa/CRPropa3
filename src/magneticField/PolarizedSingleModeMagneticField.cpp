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
	if (e_1.getR() > 0) {
		unitVector_1 = e_1 / e_1.getR();
	}
	else {
		throw std::runtime_error("Vector e_1 cannot be zero.");
	}

	Vector3d unitVector_2;
	if (e_2.getR() > 0) {
		unitVector_2 = e_2 / e_2.getR();
	}
	else {
                throw std::runtime_error("Vector e_2 cannot be zero.");
        }

	Vector3d wavevector = e_2.cross(e_1);
	if (wavevector.getR() > 0) {
                wavevector = wavevector / wavevector.getR();
        }
        else {
                throw std::runtime_error("e_1 cannot be parallel to e_2.");
        }

	if (wavelength > 0) {
		wavevector = wavevector/wavevector.getR() * 2 * M_PI / wavelength;
	}
        else {
                throw std::runtime_error("The correlation length cannot be zero.");
        }

	return B_0 * (unitVector_1 * cos(wavevector.dot(deltar)) + unitVector_2 * sigma * sin(wavevector.dot(deltar)));
}

Vector3d PolarizedSingleModeMagneticField::getBrmsMaxHelFracField(const Vector3d &pos, const double &fH) const {
	if (abs(fH > 1)) {
		throw std::runtime_error("The value of sigma has to be in the range [-1;1].");
	}
	if (fH == 0) {
		getField(pos, 0);
	}
	else {
		return sqrt(1 + sqrt(1 - fH * fH))*getField(pos,(1 - sqrt(1 - fH * fH))/fH);
	}
}

Vector3d PolarizedSingleModeMagneticField::getOrthogonalBrmsMaxHelFracField(const Vector3d &pos, const double &fH) const {
	if (e_1.dot(e_2) != 0) {
                throw std::runtime_error("e_1 and e_2 have to be orthogonal to each other.");
        }
        else {
                return getBrmsMaxHelFracField(pos, fH);
        }	
}

Vector3d PolarizedSingleModeMagneticField::getGeneralOrthogonalEllipticSingleModeB0MagneticField(const Vector3d &pos, const double &sigma) const {
        if (e_1.dot(e_2) != 0) {
                throw std::runtime_error("e_1 and e_2 have to be orthogonal to each other.");
        }
        else {
                return getField(pos, sigma);
        }
}

Vector3d PolarizedSingleModeMagneticField::getGeneralOrthogonalEllipticSingleModeBrmsMagneticField(const Vector3d &pos, const double &sigma) const {
    	return sqrt(2 / (1 + sigma * sigma)) * getGeneralOrthogonalEllipticSingleModeB0MagneticField(pos, sigma);
}

Vector3d PolarizedSingleModeMagneticField::getSpecialOrthogonalEllipticSingleModeB0MagneticField(const Vector3d &pos, const double &sigma) const {
	if (abs(sigma) > 1) {
		throw std::runtime_error("The value of sigma has to be in the range [-1;1].");
	}
	else {
		return getGeneralOrthogonalEllipticSingleModeB0MagneticField(pos, sigma);
	}
}

Vector3d PolarizedSingleModeMagneticField::getSpecialOrthogonalEllipticSingleModeBrmsMagneticField(const Vector3d &pos, const double &sigma) const {
	return sqrt(2 / (1 + sigma * sigma)) * getSpecialOrthogonalEllipticSingleModeB0MagneticField(pos, sigma);
}

Vector3d PolarizedSingleModeMagneticField::getCircularSingleModeMagneticField(const Vector3d &pos, const double &sigma) const {
        if (sigma == 0) {
                throw std::runtime_error("The value of sigma cannot be zero.");
        }
        else {
                return getGeneralOrthogonalEllipticSingleModeB0MagneticField(pos, sigma/abs(sigma));
        }
}

Vector3d PolarizedSingleModeMagneticField::getLinearSingleModeB0MagneticField(const Vector3d &pos) const {
	return getGeneralOrthogonalEllipticSingleModeB0MagneticField(pos, 0);
}

Vector3d PolarizedSingleModeMagneticField::getLinearSingleModeBrmsMagneticField(const Vector3d &pos) const {
	return sqrt(2) * getLinearSingleModeB0MagneticField(pos);
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

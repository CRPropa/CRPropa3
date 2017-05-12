#include "crpropa/magneticField/ArchmedeanSpiralField.h"

namespace crpropa {

ArchmedeanSpiralField::ArchmedeanSpiralField(double B_0, double R_0, double Omega, double V_w) {
	setB0(B_0);
	setR0(R_0);
	setOmega(Omega);
	setVw(V_w);
}

Vector3d ArchmedeanSpiralField::getField(const Vector3d &pos) const {
	
	double r = pos.getR();
	double theta = pos.getTheta();
	double phi =pos.getPhi();
	
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);

	Vector3d B(0.);

	// radial direction
	double C1 = R_0*R_0/r/r;
	B.x += C1 * cos_phi * sin_theta;
	B.y += C1 * sin_phi * sin_theta;
	B.z += C1 * cos_theta;
	
	// azimuthal direction	
	double C2 = - (Omega*R_0*R_0*sin_theta) / (r*V_w);
	B.x += C2 * (-sin_phi) * sin_theta;
	B.y += C2 * cos_phi * sin_theta;

	// magnetic field switch at z = 0
	if (pos.z<0.) {
		B *= -1;
	}

	// overall scaling
	B *= B_0;
	

	return B;
}

void ArchmedeanSpiralField::setR0(double R) {
	R_0 = R;
	return;
}

void ArchmedeanSpiralField::setB0(double B) {
	B_0 = B;
	return;
}

void ArchmedeanSpiralField::setOmega(double Om) {
	Omega = Om;
	return;
}

void ArchmedeanSpiralField::setVw(double v) {
	V_w = v;
	return;
}


double ArchmedeanSpiralField::getR0() const {
	return R_0;
}

double ArchmedeanSpiralField::getB0() const{
	return B_0;
}

double ArchmedeanSpiralField::getOmega() const{
	return Omega;
}

double ArchmedeanSpiralField::getVw() const {
	return V_w;
}


} //end namespace crpropa

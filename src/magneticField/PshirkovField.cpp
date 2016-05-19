#include "crpropa/magneticField/PshirkovField.h"
#include "crpropa/Units.h"

namespace crpropa {

PshirkovField::PshirkovField() {
	setUseBSS(true);
	setUseHalo(true);

	// disk parameters
	d = - 0.6 * kpc;
	R_sun = 8.5 * kpc;
	R_c = 5.0 * kpc;
	z0 = 1.0 * kpc;
	B0 = 2.0 * muG;

	// halo parameters
	z0_H = 1.3 * kpc;
	R0_H = 8.0 * kpc;
	B0_H = 4.0 * muG;
	B0_H_ASS = 2.0 * muG;
	B0_H_BSS = 4.0 * muG;
	z11_H = 0.25 * kpc;
	z12_H = 0.4 * kpc;
}

void PshirkovField::setUseASS(bool use) {
	useASS = use;
	if ((use) and (useBSS)) {
		std::cout << "PshirkovField: Disk field changed to ASS" << std::endl;
		useBSS = false;
	}

	pitch = -5.0 * M_PI / 180;
	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
	theta = cos_pitch / sin_pitch * log(1 + d / R_sun) - M_PI / 2;
	cos_theta = cos(theta);
}

void PshirkovField::setUseBSS(bool use) {
	useBSS = use;
	if ((use) and (useASS)) {
		std::cout << "PshirkovField: Disk field changed to BSS" << std::endl;
		useASS = false;
	}

	pitch = -6.0 * M_PI / 180;
	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
	theta = cos_pitch / sin_pitch * log(1 + d / R_sun) - M_PI / 2;
	cos_theta = cos(theta);
}

void PshirkovField::setUseHalo(bool use) {
	useHalo = use;
}

bool PshirkovField::isUsingASS() {
	return useASS;
}

bool PshirkovField::isUsingBSS() {
	return useBSS;
}

bool PshirkovField::isUsingHalo() {
	return useHalo;
}

Vector3d PshirkovField::getDiskField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (pos.getR() < 0.1 * kpc)
		return b; // no disk field for distances to galactic center < 0.1 kpc

	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
	double phi = pos.getPhi(); // azimuth
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);
	b.x += sin_pitch * cos_phi - cos_pitch * sin_phi;
	b.y += sin_pitch * sin_phi + cos_pitch * cos_phi;

	double bMag = cos(phi - cos_pitch / sin_pitch * log(r / R_sun) + theta);
	if (useASS)
		bMag = fabs(bMag);
	bMag *= exp(-fabs(pos.z) / z0);

	if (r <= R_c)
		bMag *= B0 * R_sun / (cos_theta * R_c);  // inner disk field
	else
		bMag *= B0 * R_sun / (cos_theta * r);  // outer disk field

	return bMag * b;
}

Vector3d PshirkovField::getHaloField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);
	double bMag = pow(1 + pow((fabs(pos.z) - z0_H) / z11_H, 2.), -1.) * r / R0_H * exp(1 - r / R0_H);

	if (pos.z > 0)
		bMag *= B0_H;  // northern halo
	else
		bMag *= - (useASS? B0_H_ASS : B0_H_BSS);  // southern halo

	double phi = pos.getPhi();
	return Vector3d(-sin(phi), cos(phi), 0) * bMag;
}

Vector3d PshirkovField::getField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (pos.getR() > 20 * kpc)
		return b; // no field for distances to galactic center > 20 kpc

	if (useHalo)
		b += getHaloField(pos);
	if ((useASS) or (useBSS))
		b += getDiskField(pos);
	return b;
}

} // namespace crpropa

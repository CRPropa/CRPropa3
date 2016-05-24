#include "crpropa/magneticField/PshirkovField.h"
#include "crpropa/Units.h"

namespace crpropa {

PshirkovField::PshirkovField() {
	// disk parameters
	d = - 0.6 * kpc;
	R_sun = 8.5 * kpc;
	R_c = 5.0 * kpc;
	z0 = 1.0 * kpc;
	B0 = 2.0 * muG;

	// halo parameters
	z0_H = 1.3 * kpc;
	R0_H = 8.0 * kpc;
	B0_Hn = 4.0 * muG;
	B0_Hs = 4.0 * muG;
	z11_H = 0.25 * kpc;
	z12_H = 0.4 * kpc;

	setUseBSS(true);
	setUseHalo(true);
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
	B0_Hs = 2.0 * muG;
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
	B0_Hs = 4.0 * muG;
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
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
	double phi = pos.getPhi(); // azimuth
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);

	Vector3d b(0.);
	b.x = sin_pitch * cos_phi - cos_pitch * sin_phi;
	b.y = sin_pitch * sin_phi + cos_pitch * cos_phi;

	double bMag = cos(phi - cos_pitch / sin_pitch * log(r / R_sun) + theta);
	if (useASS)
		bMag = fabs(bMag);
	bMag *= B0 * R_sun / cos_theta * exp(-fabs(pos.z) / z0);
	bMag /= (r <= R_c ? R_c : r);

	return bMag * b;
}

Vector3d PshirkovField::getHaloField(const Vector3d& pos) const {
	double bMag = (pos.z > 0 ? B0_Hn : - B0_Hs);
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);
	bMag *= r / R0_H * exp(1 - r / R0_H) / (1 + pow((fabs(pos.z) - z0_H) / z11_H, 2.));
	double phi = pos.getPhi();
	return bMag * Vector3d(-sin(phi), cos(phi), 0);
}

Vector3d PshirkovField::getField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (useHalo)
		b += getHaloField(pos);
	if ((useASS) or (useBSS))
		b += getDiskField(pos);
	return b;
}

} // namespace crpropa

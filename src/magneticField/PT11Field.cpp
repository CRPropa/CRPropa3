#include "crpropa/magneticField/PT11Field.h"
#include "crpropa/Units.h"

#include <algorithm>

namespace crpropa {

PT11Field::PT11Field() : useASS(false), useBSS(true), useHalo(true) {
	// disk parameters
	d = - 0.6 * kpc;
	R_sun = 8.5 * kpc;
	R_c = 5.0 * kpc;
	z0_D = 1.0 * kpc;
	B0_D = 2.0 * muG;

	// halo parameters
	z0_H = 1.3 * kpc;
	R0_H = 8.0 * kpc;
	B0_Hn = 4.0 * muG;
	B0_Hs = 4.0 * muG;
	z11_H = 0.25 * kpc;
	z12_H = 0.4 * kpc;

	// set BSS specific parameters
	setUseBSS(true);
}

void PT11Field::setUseASS(bool use) {
	useASS = use;
	if (not(use))
		return;

	if (useBSS) {
		std::cout << "PT11Field: Disk field changed to ASS" << std::endl;
		useBSS = false;
	}

	pitch = -5.0 * M_PI / 180;
	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
	theta = cos_pitch / sin_pitch * log(1 + d / R_sun) - M_PI / 2;
	cos_theta = cos(theta);
	B0_Hs = 2.0 * muG;
}

void PT11Field::setUseBSS(bool use) {
	useBSS = use;
	if (not(use))
		return;

	if (useASS) {
		std::cout << "PT11Field: Disk field changed to BSS" << std::endl;
		useASS = false;
	}

	pitch = -6.0 * M_PI / 180;
	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
	theta = cos_pitch / sin_pitch * log(1 + d / R_sun) - M_PI / 2;
	cos_theta = cos(theta);
	B0_Hs = 4.0 * muG;
}

void PT11Field::setUseHalo(bool use) {
	useHalo = use;
}

bool PT11Field::isUsingASS() {
	return useASS;
}

bool PT11Field::isUsingBSS() {
	return useBSS;
}

bool PT11Field::isUsingHalo() {
	return useHalo;
}

Vector3d PT11Field::getField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius

	double cos_phi = pos.x / r;
	double sin_phi = pos.y / r;

	Vector3d b(0.);

	// disk field
	if ((useASS) or (useBSS)) {
		// PT11 paper has B * cos(p) but this seems because they define azimuth clockwise, while we have anticlockwise.
		// see Tinyakov 2002 APh 18,165: "local field points to l=90+p" so p=-5 deg gives l=85 and hence clockwise from above.
		// so to get local B clockwise in our system, need minus (like Sun etal).
		// Ps base their system on Han and Qiao 1994 A&A 288,759 which has a diagram with azimuth clockwise, hence confirmed.
		double phi = -pos.getPhi();  // azimuth; since PT11 paper uses opposite convention for phi
		double bMag = cos(phi - cos_pitch / sin_pitch * log(r / R_sun) + theta);
		double cos_pitch_pt = - cos_pitch; // azimuthal field in direction of increasing azimuth angle theta
		b.x = sin_pitch * cos_phi - cos_pitch_pt * sin_phi;
		b.y = sin_pitch * sin_phi + cos_pitch_pt * cos_phi;
		if (useASS)
			bMag = fabs(bMag);
		bMag *= B0_D * R_sun / std::max(r, R_c) / cos_theta * exp(-fabs(pos.z) / z0_D);
		b *= bMag;
	}

	// halo field
	if (useHalo) {
		double bMag = (pos.z > 0 ? B0_Hn : - B0_Hs);
		double z1 = (fabs(pos.z) < z0_H ? z11_H : z12_H);
		bMag *= r / R0_H * exp(1 - r / R0_H) / (1 + pow((fabs(pos.z) - z0_H) / z1, 2.));
		b += bMag * Vector3d(- sin_phi, cos_phi, 0);
	}

	return b;
}

} // namespace crpropa

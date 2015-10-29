#include "crpropa/magneticField/PshirkovField.h"
#include "crpropa/Units.h"

#include <iostream>

namespace crpropa {

PshirkovField::PshirkovField() {
	useBSS = true;
	useBHMHalo = true;

	// disk parameters
	pitch_ASS = -5.0 * M_PI / 180;
	pitch_BSS = -6.0 * M_PI / 180;
	b_ASS = 1. / tan(pitch_ASS);
	b_BSS = 1. / tan(pitch_BSS);
	d = - 0.6 * kpc;
	R_sun = 8.5 * kpc;
	R_c = 5.0 * kpc;
	theta_ASS = b_ASS * log(1+d/R_sun) - M_PI / 2;
	theta_BSS = b_BSS * log(1+d/R_sun) - M_PI / 2;
	z0 = 1.0 * kpc;
	B0 = 2.0 * muG;

	// halo (north)
	z0_n = 1.3 * kpc;
	R0_n = 8.0 * kpc;
	B0_n = 4.0 * muG;
	z11_n = 0.25 * kpc;
	z12_n = 0.4 * kpc;

	// halo (south)
	z0_s = 1.3 * kpc;
	R0_s = 8.0 * kpc;
	B0_s_ASS = 2.0 * muG;
	B0_s_BSS = 4.0 * muG;
	z11_s = 0.25 * kpc;
	z12_s = 0.4 * kpc;
}

void PshirkovField::setUseASS(bool use) {
	useASS = use;
	if ((use) and (useBSS)) {
		std::cout << "PshirkovField: Disk field changed to ASS" << std::endl;
		useBSS = false;
	}
}

void PshirkovField::setUseBSS(bool use) {
	useBSS = use;
	if ((use) and (useASS)) {
		std::cout << "PshirkovField: Disk field changed to BSS" << std::endl;
		useASS = false;
	}
}

void PshirkovField::setUseBHMHalo(bool use) {
	useBHMHalo = use;
}

bool PshirkovField::isUsingASS() {
	return useASS;
}

bool PshirkovField::isUsingBSS() {
	return useBSS;
}

bool PshirkovField::isUsingBHMHalo() {
	return useBHMHalo;
}

double PshirkovField::getSouthernFieldstrength() const{
	if (useASS)
		return B0_s_ASS;
	else if (useBSS)
		return B0_s_BSS;
}


Vector3d PshirkovField::getASSField(const Vector3d& pos) const {
	Vector3d b(0.);

	double r = sqrt(pos.x * pos.x + pos.y * pos.y);                        // in-plane radius
	double R = pos.getR();                                                    // distance to galactic center
	if ((r < 0.1 * kpc) or (r > 20 * kpc)) {
		return b;                                                            // 0 field for d < 1 kpc or d > 20 kpc
	}
	double phi = pos.getPhi();                                                // azimuth
	double bMag;

	// inner disk field
	if ((r > 0.1 * kpc) and (r <= R_c)) {
		bMag = B0 * R_sun / (cos(theta_ASS) * R_c);
		double bGes = bMag * fabs(cos(phi - b_BSS * log(r / R_sun) + theta_BSS)) * exp(-fabs(pos.z) / z0);
		b.x += bGes * (sin(pitch_BSS) * cos(phi) - cos(pitch_BSS) * sin(phi));
		b.y += bGes * (sin(pitch_BSS) * sin(phi) + cos(pitch_BSS) * cos(phi));
	}
	// outer disk field
	else if ((r > R_c) and (r <= 20 * kpc)) {
		bMag = B0 * R_sun / (cos(theta_ASS) * r);
		double bGes = bMag * fabs(cos(phi - b_BSS * log(r / R_sun) + theta_BSS)) * exp(-fabs(pos.z) / z0);
		b.x += bGes * (sin(pitch_BSS) * cos(phi) - cos(pitch_BSS) * sin(phi));
		b.y += bGes * (sin(pitch_BSS) * sin(phi) + cos(pitch_BSS) * cos(phi));
	}
	return b;
}

Vector3d PshirkovField::getBSSField(const Vector3d& pos) const {
	Vector3d b(0.);

	double r = sqrt(pos.x * pos.x + pos.y * pos.y);                        // in-plane radius
	double R = pos.getR();                                                    // distance to galactic center
	if ((r <= 0.1 * kpc) or (r > 20 * kpc)) {
		return b;                                                            // 0 field for d < 1 kpc or d > 20 kpc
	}
	double phi = pos.getPhi();                                                // azimuth
	double bMag;

	// inner disk field
	if ((r > 0.1 * kpc) and (r <= R_c)) {
		bMag = B0 * R_sun / (cos(theta_BSS) * R_c);
		double bGes = bMag * cos(phi - b_BSS * log(r / R_sun) + theta_BSS) * exp(-fabs(pos.z) / z0);
		b.x += bGes * (sin(pitch_BSS) * cos(phi) - cos(pitch_BSS) * sin(phi));
		b.y += bGes * (sin(pitch_BSS) * sin(phi) + cos(pitch_BSS) * cos(phi));
	}
		// outer disk field
	else if ((r > R_c) and (r <= 20 * kpc)) {
		bMag = B0 * R_sun / (cos(theta_BSS) * r);
		double bGes = bMag * cos(phi - b_BSS * log(r / R_sun) + theta_BSS) * exp(-fabs(pos.z) / z0);
		b.x += bGes * (sin(pitch_BSS) * cos(phi) - cos(pitch_BSS) * sin(phi));
		b.y += bGes * (sin(pitch_BSS) * sin(phi) + cos(pitch_BSS) * cos(phi));
	}
	return b;
}

Vector3d PshirkovField::getBHMHaloField(const Vector3d& pos) const {
	Vector3d b(0.);
	// halo field (Basic Halo Model from Sun et al. 2008)
	double bMag_H;
	double phi = pos.getPhi();
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);
	if (fabs(pos.z) < z0_n) {
		if (pos.z >= 0) {
			bMag_H = B0_n * pow(1 + pow((fabs(pos.z) - z0_n) / z11_n, 2.), -1.) * r / R0_n * exp(1 - r / R0_n);
		}
		else {
			bMag_H = - getSouthernFieldstrength() * pow(1 + pow((fabs(pos.z) - z0_s) / z11_s, 2.), -1.) * r / R0_s * exp(1 - r / R0_s);
		}
	}
	else {
		if (pos.z >= 0) {
			bMag_H = B0_n * pow(1 + pow((fabs(pos.z) - z0_n) / z12_n, 2.), -1.) * r / R0_n * exp(1 - r / R0_n);
		}
		else {
			bMag_H = - getSouthernFieldstrength() * pow(1 + pow((fabs(pos.z) - z0_s) / z12_s, 2.), -1.) * r / R0_s * exp(1 - r / R0_s);
		}
	}
	b.x += - bMag_H * sin(phi);
	b.y += bMag_H * cos(phi);
	return b;
}

Vector3d PshirkovField::getField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (useBHMHalo)
		b += getBHMHaloField(pos);
	if (useASS)
		b += getASSField(pos);
	else if (useBSS)
		b += getBSSField(pos);
	return b;
}

} // namespace crpropa

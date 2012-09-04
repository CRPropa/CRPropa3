#include "mpc/magneticField/JF2012Field.h"
#include "mpc/Units.h"
#include "mpc/GridTools.h"

namespace mpc {

double logisticFunction(double x, double x0, double w) {
	return 1. / (1. + exp(-2. * (fabs(x) - x0) / w));
}

JF2012Field::JF2012Field() {
	pitch = 11.5 * M_PI / 180;
	sinPitch = sin(pitch);
	cosPitch = cos(pitch);
	tan90MinusPitch = tan(M_PI / 2 - pitch);

	bRing = 0.1 * muG;
	hDisk = 0.40 * kpc;
	wDisk = 0.27 * kpc;

	bDisk[0] = 0.1 * muG;
	bDisk[1] = 3.0 * muG;
	bDisk[2] = -0.9 * muG;
	bDisk[3] = -0.8 * muG;
	bDisk[4] = -2.0 * muG;
	bDisk[5] = -4.2 * muG;
	bDisk[6] = 0.0 * muG;
	bDisk[7] = 2.7 * muG;

	rArms[0] = 5.1 * kpc;
	rArms[1] = 6.3 * kpc;
	rArms[2] = 7.1 * kpc;
	rArms[3] = 8.3 * kpc;
	rArms[4] = 9.8 * kpc;
	rArms[5] = 11.4 * kpc;
	rArms[6] = 12.7 * kpc;
	rArms[7] = 15.5 * kpc;

	bNorth = 1.4 * muG;
	bSouth = -1.1 * muG;
	rNorth = 9.22 * kpc;
	rSouth = 17 * kpc;
	wHalo = 0.20 * kpc;
	z0 = 5.3 * kpc;

	bX = 4.6 * muG;
	thetaX0 = 49.0 * M_PI / 180;
	sinThetaX0 = sin(thetaX0);
	cosThetaX0 = cos(thetaX0);
	tanThetaX0 = tan(thetaX0);
	rXc = 4.8 * kpc;
	rX = 2.9 * kpc;
}

Vector3d JF2012Field::getField(const Vector3d& pos) const {
	Vector3d b(0.);

	if (pos.getMag() < 1 * kpc)
		return b; // no field if distance to GC < 1 kpc

	double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
	if (r > 20 * kpc)
		return b; // no field if in-plane distance to GC > 20 kpc

	double phi = pos.getPhi(); // azimuth
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);

	double lfDisk = logisticFunction(pos.z, hDisk, wDisk);

	// disk field
	if (r > 3 * kpc) {
		if (r < 5 * kpc) {
			// molecular ring
			double bMagD = bRing * (5 * kpc / r) * (1 - lfDisk);
			b.x += -bMagD * sinPhi;
			b.y += bMagD * cosPhi;

		} else if (r < 20 * kpc) {
			// spiral region
			double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

			double bMagD;
			for (int i = 7; i >= 0; i--)
				if (r_negx < rArms[i])
					bMagD = bDisk[i];

			bMagD *= (5 * kpc / r) * (1 - lfDisk);
			b.x += bMagD * (sinPitch * cosPhi - cosPitch * sinPhi);
			b.y += bMagD * (sinPitch * sinPhi + cosPitch * cosPhi);
		}
	}

	// toroidal halo field
	double bMagH = exp(-fabs(pos.z) / z0) * lfDisk;
	if (pos.z >= 0)
		bMagH *= bNorth * (1 - logisticFunction(r, rNorth, wHalo));
	else
		bMagH *= bSouth * (1 - logisticFunction(r, rSouth, wHalo));
	b.x += -bMagH * sinPhi;
	b.y += bMagH * cosPhi;

	// X-halo field
	double bMagX;
	double sinThetaX, cosThetaX;
	double rp;
	double rc = rXc + fabs(pos.z) / tanThetaX0;
	if (r < rc) {
		// varying elevation region
		rp = r * rXc / rc;
		bMagX = bX * exp(-rp / rX) * pow(rp / r, 2.);
		double thetaX = atan2(fabs(pos.z), (r - rp));
		sinThetaX = sin(thetaX);
		cosThetaX = cos(thetaX);
	} else {
		// constant elevation region
		rp = r - fabs(pos.z) / tanThetaX0;
		bMagX = bX * exp(-rp / rX) * (rp / r);
		sinThetaX = sinThetaX0;
		cosThetaX = cosThetaX0;
	}
	double zsign = pos.z < 0 ? -1 : 1;
	b.x += zsign * bMagX * cosThetaX * cosPhi;
	b.y += zsign * bMagX * cosThetaX * sinPhi;
	b.z += bMagX * sinThetaX;

	return b;
}

JF2012TurbulentField::JF2012TurbulentField() {
	// correlation length = 0.060 * kpc;
	r0 = 10.65363798 * kpc;
	z0 = 3.89915783 * kpc;

	pitch = 11.5 * M_PI / 180;
	tan90MinusPitch = tan(M_PI / 2 - pitch);

	b5 = 5.35262527 * muG;
	z0S = 0.58119582 * kpc;

	bDisk[0] = 7.94373074 * muG;
	bDisk[1] = 4.4140519 * muG;
	bDisk[2] = 6.72729298 * muG;
	bDisk[3] = 4.69436688 * muG;
	bDisk[4] = 1.35254503 * muG;
	bDisk[5] = 10.91095824;
	bDisk[6] = 24.74973782;
	bDisk[7] = 6.86975446;

	rArms[0] = 5.1 * kpc;
	rArms[1] = 6.3 * kpc;
	rArms[2] = 7.1 * kpc;
	rArms[3] = 8.3 * kpc;
	rArms[4] = 9.8 * kpc;
	rArms[5] = 11.4 * kpc;
	rArms[6] = 12.7 * kpc;
	rArms[7] = 15.5 * kpc;
}

double JF2012TurbulentField::getModulation(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);
	if (r > 20 * kpc)
		return 0; // no field if in-plane distance to GC > 20 kpc

	double phi = pos.getPhi(); // azimuth

	// disk field
	double spiralNorm;
	if (r < 5 * kpc) {
		spiralNorm = b5;
	} else if (r < 20 * kpc) {
		double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

		for (int i = 7; i >= 0; i--)
			if (r_negx < rArms[i])
				spiralNorm = bDisk[i];

		spiralNorm *= (5 * kpc / r);
	}
	spiralNorm *= exp(-0.5 * pow(pos.z / z0S, 2));

	double smoothNorm = 1;
	smoothNorm *= exp(-1. * fabs(pos.z) / z0);
	smoothNorm *= exp(-1. * fabs(r) / r0);

	return sqrt(smoothNorm * smoothNorm + spiralNorm * spiralNorm);
}

Vector3d JF2012TurbulentField::getField(const Vector3d& pos) const {
	Vector3d b = grid->interpolate(pos);
	return b * getModulation(pos);
}

}

#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

namespace crpropa {

double logisticFunction(double x, double x0, double w) {
	return 1. / (1. + exp(-2. * (fabs(x) - x0) / w));
}

JF12Field::JF12Field() {
	useRegular = true;
	useStriated = false;
	useTurbulent = false;

	// spiral arm parameters
	pitch = 11.5 * M_PI / 180;
	sinPitch = sin(pitch);
	cosPitch = cos(pitch);
	tan90MinusPitch = tan(M_PI / 2 - pitch);

	rArms[0] = 5.1 * kpc;
	rArms[1] = 6.3 * kpc;
	rArms[2] = 7.1 * kpc;
	rArms[3] = 8.3 * kpc;
	rArms[4] = 9.8 * kpc;
	rArms[5] = 11.4 * kpc;
	rArms[6] = 12.7 * kpc;
	rArms[7] = 15.5 * kpc;

	// regular field parameters
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

	// striated field parameter
	sqrtbeta = sqrt(1.36);

	// turbulent field parameters
	bDiskTurb[0] = 10.81 * muG;
	bDiskTurb[1] = 6.96 * muG;
	bDiskTurb[2] = 9.59 * muG;
	bDiskTurb[3] = 6.96 * muG;
	bDiskTurb[4] = 1.96 * muG;
	bDiskTurb[5] = 16.34 * muG;
	bDiskTurb[6] = 37.29 * muG;
	bDiskTurb[7] = 10.35 * muG;

	bDiskTurb5 = 7.63 * muG;
	zDiskTurb = 0.61 * kpc;

	bHaloTurb = 4.68 * muG;
	rHaloTurb = 10.97 * kpc;
	zHaloTurb = 2.84 * kpc;
}

void JF12Field::randomStriated(int seed) {
	useStriated = true;
	int N = 100;
	striatedGrid = new ScalarGrid(Vector3d(0.), N, 0.1 * kpc);

	Random random;
	if (seed != 0)
		random.seed(seed);

	for (int ix = 0; ix < N; ix++)
		for (int iy = 0; iy < N; iy++)
			for (int iz = 0; iz < N; iz++) {
				float &f = striatedGrid->get(ix, iy, iz);
				f = round(random.rand()) * 2 - 1;
			}
}

#ifdef CRPROPA_HAVE_FFTW3F
void JF12Field::randomTurbulent(int seed) {
	useTurbulent = true;
	// turbulent field with Kolmogorov spectrum, B_rms = 1 and Lc = 60 parsec
	turbulentGrid = new VectorGrid(Vector3d(0.), 256, 4 * parsec);
	initTurbulence(turbulentGrid, 1, 8 * parsec, 272 * parsec, -11./3., seed);
}
#endif

void JF12Field::setStriatedGrid(ref_ptr<ScalarGrid> grid) {
	useStriated = true;
	striatedGrid = grid;
}

void JF12Field::setTurbulentGrid(ref_ptr<VectorGrid> grid) {
	useTurbulent = true;
	turbulentGrid = grid;
}

ref_ptr<ScalarGrid> JF12Field::getStriatedGrid() {
	return striatedGrid;
}

ref_ptr<VectorGrid> JF12Field::getTurbulentGrid() {
	return turbulentGrid;
}

void JF12Field::setUseRegular(bool use) {
	useRegular = use;
}

void JF12Field::setUseStriated(bool use) {
	if ((use) and (striatedGrid)) {
		std::cout << "JF12Field: No striated field set: ignored" << std::endl;
		return;
	}
	useStriated = use;
}

void JF12Field::setUseTurbulent(bool use) {
	if ((use) and (turbulentGrid)) {
		std::cout << "JF12Field: No turbulent field set: ignored" << std::endl;
		return;
	}
	useTurbulent = use;
}

bool JF12Field::isUsingRegular() {
	return useRegular;
}

bool JF12Field::isUsingStriated() {
	return useStriated;
}

bool JF12Field::isUsingTurbulent() {
	return useTurbulent;
}

Vector3d JF12Field::getRegularField(const Vector3d& pos) const {
	Vector3d b(0.);

	double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
	double d = pos.getR(); // distance to galactic center
	if ((d < 1 * kpc) or (d > 20 * kpc))
		return b; // 0 field for d < 1 kpc or d > 20 kpc

	double phi = pos.getPhi(); // azimuth
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);

	double lfDisk = logisticFunction(pos.z, hDisk, wDisk);

	// disk field
	if (r > 3 * kpc) {
		double bMag;
		if (r < 5 * kpc) {
			// molecular ring
			bMag = bRing * (5 * kpc / r) * (1 - lfDisk);
			b.x += -bMag * sinPhi;
			b.y += bMag * cosPhi;

		} else {
			// spiral region
			double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

			for (int i = 7; i >= 0; i--)
				if (r_negx < rArms[i])
					bMag = bDisk[i];

			bMag *= (5 * kpc / r) * (1 - lfDisk);
			b.x += bMag * (sinPitch * cosPhi - cosPitch * sinPhi);
			b.y += bMag * (sinPitch * sinPhi + cosPitch * cosPhi);
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

	// poloidal halo field
	double bMagX;
	double sinThetaX, cosThetaX;
	double rp;
	double rc = rXc + fabs(pos.z) / tanThetaX0;
	if (r < rc) {
		// varying elevation region
		rp = r * rXc / rc;
		bMagX = bX * exp(-1 * rp / rX) * pow(rXc / rc, 2.);
		double thetaX = atan2(fabs(pos.z), (r - rp));
		if (pos.z == 0)
			thetaX = M_PI / 2.;
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

Vector3d JF12Field::getStriatedField(const Vector3d& pos) const {
	return (getRegularField(pos)
			* (1. + sqrtbeta * striatedGrid->closestValue(pos)));
}

double JF12Field::getTurbulentStrength(const Vector3d& pos) const {
	if (pos.getR() > 20 * kpc)
		return 0;

	double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
	double phi = pos.getPhi(); // azimuth

	// disk
	double bDisk = 0;
	if (r < 5 * kpc) {
		bDisk = bDiskTurb5;
	} else {
		// spiral region
		double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
		if (r_negx > rArms[7])
			r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

		for (int i = 7; i >= 0; i--)
			if (r_negx < rArms[i])
				bDisk = bDiskTurb[i];

		bDisk *= (5 * kpc) / r;
	}
	bDisk *= exp(-0.5 * pow(pos.z / zDiskTurb, 2));

	// halo
	double bHalo = bHaloTurb * exp(-r / rHaloTurb)
			* exp(-0.5 * pow(pos.z / zHaloTurb, 2));

	// modulate turbulent field
	return sqrt(pow(bDisk, 2) + pow(bHalo, 2));
}

Vector3d JF12Field::getTurbulentField(const Vector3d& pos) const {
	return (turbulentGrid->interpolate(pos) * getTurbulentStrength(pos));
}

Vector3d JF12Field::getField(const Vector3d& pos) const {
	Vector3d b(0.);
	if (useTurbulent)
		b += getTurbulentField(pos);
	if (useStriated)
		b += getStriatedField(pos);
	else if (useRegular)
		b += getRegularField(pos);
	return b;
}

} // namespace crpropa

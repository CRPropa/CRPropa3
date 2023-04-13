#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"
#include "crpropa/Random.h"

namespace crpropa {

JF12Field::JF12Field() {
	useRegularField = true;
	useStriatedField = false;
	useTurbulentField = false;
	useDiskField = true;
	useToroidalHaloField = true;
	useXField = true;

	// spiral arm parameters
	pitch = 11.5 * M_PI / 180;
	sinPitch = sin(pitch);
	cosPitch = cos(pitch);
	tanPitch = tan(pitch);
	cotPitch =  1. / tanPitch;
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
	cotThetaX0 = 1. / tanThetaX0;
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
	useStriatedField = true;
	int N = 100;
	striatedGrid = new Grid1f(Vector3d(0.), N, 0.1 * kpc);

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
	useTurbulentField = true;
	// turbulent field with Kolmogorov spectrum, B_rms = 1 (will be scaled) and Lc = 60 parsec, and 256 grid points.
	// Note that the inertial range of the turbulence is less than 2 orders of magnitude.
    const double lMin = 8 * parsec;
    const double lMax = 272 * parsec;
    const double Brms = 1;
    const double spacing = 4 * parsec;
    const double grid_n = 256;

    auto spectrum = SimpleTurbulenceSpectrum(Brms, lMin, lMax);
    auto gp = GridProperties(Vector3d(0.), grid_n, spacing);
    auto tf = SimpleGridTurbulence(spectrum, gp, seed);
    turbulentGrid = tf.getGrid();

}
#endif

void JF12Field::setStriatedGrid(ref_ptr<Grid1f> grid) {
	useStriatedField = true;
	striatedGrid = grid;
}

void JF12Field::setTurbulentGrid(ref_ptr<Grid3f> grid) {
	useTurbulentField = true;
	turbulentGrid = grid;
}

ref_ptr<Grid1f> JF12Field::getStriatedGrid() {
	return striatedGrid;
}

ref_ptr<Grid3f> JF12Field::getTurbulentGrid() {
	return turbulentGrid;
}

void JF12Field::setUseRegularField(bool use) {
	useRegularField = use;
}

void JF12Field::setUseDiskField(bool use) {
	useDiskField = use;
}

void JF12Field::setUseToroidalHaloField(bool use) {
	useToroidalHaloField = use;
}

void JF12Field::setUseXField(bool use) {
	useXField = use;
}

void JF12Field::setUseStriatedField(bool use) {
	if ((use) and !(striatedGrid)) {
		KISS_LOG_WARNING << "JF12Field: No striated field set: ignored. Run e.g. randomStriated().";
		return;
	}
	useStriatedField = use;
}

void JF12Field::setUseTurbulentField(bool use) {
	if ((use) and !(turbulentGrid)) {
		KISS_LOG_WARNING << "JF12Field: No turbulent field set: ignored. Run e.g. randomTurbulent().";
		return;
	}
	useTurbulentField = use;
}

bool JF12Field::isUsingRegularField() {
	return useRegularField;
}

bool JF12Field::isUsingDiskField() {
	return useDiskField;
}

bool JF12Field::isUsingToroidalHaloField() {
	return useToroidalHaloField;
}

bool JF12Field::isUsingXField() {
	return useXField;
}

bool JF12Field::isUsingStriatedField() {
	return useStriatedField;
}

bool JF12Field::isUsingTurbulentField() {
	return useTurbulentField;
}

double JF12Field::logisticFunction(const double& x, const double& x0, const double& w) const {
	return 1. / (1. + exp(-2. * (fabs(x) - x0) / w));
}

Vector3d JF12Field::getRegularField(const Vector3d& pos) const {
	Vector3d b(0.);

	double d = pos.getR(); // distance to galactic center

	if (d < 20 * kpc) {
		double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
		double phi = pos.getPhi(); // azimuth
		double sinPhi = sin(phi);
		double cosPhi = cos(phi);

		b += getDiskField(r, pos.z, phi, sinPhi, cosPhi);
		b += getToroidalHaloField(r, pos.z, sinPhi, cosPhi);
		b += getXField(r, pos.z, sinPhi, cosPhi);
	}

	return b;
}

Vector3d JF12Field::getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const {
	Vector3d b(0.);
	if (useDiskField) {
		double lfDisk = logisticFunction(z, hDisk, wDisk);
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
	}
	return b;
}

Vector3d JF12Field::getToroidalHaloField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const {
	Vector3d b(0.);

	if (useToroidalHaloField && (r * r + z * z > 1 * kpc * kpc)){

		double lfDisk = logisticFunction(z, hDisk, wDisk);
		double bMagH = exp(-fabs(z) / z0) * lfDisk;

		if (z >= 0)
			bMagH *= bNorth * (1 - logisticFunction(r, rNorth, wHalo));
		else
			bMagH *= bSouth * (1 - logisticFunction(r, rSouth, wHalo));
		b.x += -bMagH * sinPhi;
		b.y += bMagH * cosPhi;
	}
	return b;
}

Vector3d JF12Field::getXField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const {
	Vector3d b(0.);

	if (useXField && (r * r + z * z > 1 * kpc * kpc)){
		double bMagX;
		double sinThetaX, cosThetaX;
		double rp;
		double rc = rXc + fabs(z) / tanThetaX0;
		if (r < rc) {
			// varying elevation region
			rp = r * rXc / rc;
			bMagX = bX * exp(-1 * rp / rX) * pow(rXc / rc, 2.);
			double thetaX = atan2(fabs(z), (r - rp));
			if (z == 0)
				thetaX = M_PI / 2.;
			sinThetaX = sin(thetaX);
			cosThetaX = cos(thetaX);
		} else {
			// constant elevation region
			rp = r - fabs(z) / tanThetaX0;
			bMagX = bX * exp(-rp / rX) * (rp / r);
			sinThetaX = sinThetaX0;
			cosThetaX = cosThetaX0;
		}
		double zsign = z < 0 ? -1 : 1;
		b.x += zsign * bMagX * cosThetaX * cosPhi;
		b.y += zsign * bMagX * cosThetaX * sinPhi;
		b.z += bMagX * sinThetaX;
	}
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
	if (useTurbulentField)
		b += getTurbulentField(pos);
	if (useStriatedField)
		b += getStriatedField(pos);
	else if (useRegularField)
		b += getRegularField(pos);
	return b;
}



PlanckJF12bField::PlanckJF12bField() : JF12Field::JF12Field(){
	// regular field parameters
	bDisk[5] = -3.5 * muG;
	bX = 1.8 * muG;

	// turbulent field parameters;
	bDiskTurb[0] = 3.12 * muG;
	bDiskTurb[1] = 6.24 * muG;
	bDiskTurb[2] = 3.12 * muG;
	bDiskTurb[3] = 6.24 * muG;
	bDiskTurb[4] = 3.12 * muG;
	bDiskTurb[5] = 6.24 * muG;
	bDiskTurb[6] = 3.12 * muG;
	bDiskTurb[7] = 6.24 * muG;

	bDiskTurb5 = 3.90 * muG;

	bHaloTurb = 7.332 * muG;
}

} // namespace crpropa

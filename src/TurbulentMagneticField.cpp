#include "mpc/magneticfield/TurbulentMagneticField.h"
#include "mpc/MersenneTwister.h"
#include "fftw3.h"

namespace mpc {

TurbulentMagneticField::TurbulentMagneticField(Vector3 origin, size_t n,
		double spacing, double Brms, double spectralIndex, double lMin,
		double lMax) {
	this->origin = origin;
	this->n = n;
	this->spacing = spacing;
	this->Brms = Brms;
	this->spectralIndex = spectralIndex;
	this->lMin = lMin / spacing;
	this->lMax = lMax / spacing;
}

TurbulentMagneticField::~TurbulentMagneticField() {
}

void TurbulentMagneticField::setSeed(unsigned int seed) {
	this->seed = seed;
}

void TurbulentMagneticField::initialize() {
	// random number generator
	MTRand mtrand;
	mtrand.seed(seed);

	// arrays to hold the complex vector components of B-field
	fftw_complex *Bx, *By, *Bz;
	Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
	By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
	Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);

	// N discrete possible wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// create field in configuration space
	int i;
	double k, theta, phase, cosPhase, sinPhase;
	double kMin = 1. / lMax;
	double kMax = 1. / lMin;
	Vector3 e1, e2, ek, b;
	Vector3 n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (int ix = 0; ix < n; ix++) {
		for (int iy = 0; iy < n; iy++) {
			for (int iz = 0; iz < n; iz++) {

				i = ix * n * n + iy * n + iz;
				ek.set(K[ix], K[iy], K[iz]);
				k = ek.mag();

				if ((k < kMin) || (k > kMax))
					continue; // wave outside of turbulent range -> B(k) = 0

				// construct an orthogonal base e1, e2, ek
				if ((ix == iy) && (iy == iz)) {
					// ek parallel to (1,1,1)
					e1.set(-1., 1., 0);
					e2.set(1., 1., -2.);
				} else {
					// ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.mag();
				e2 /= e2.mag();

				// random orientation perpendicular to k
				theta = 2 * M_PI * mtrand.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// gaussian amplitude weighted with k^alpha/2
				b *= mtrand.randNorm() * pow(k, spectralIndex / 2.);

				// uniform random phase
				phase = 2 * M_PI * mtrand.rand();
				cosPhase = cos(phase); // real part
				sinPhase = sin(phase); // imaginary part

				Bx[i][0] = b.x() * cosPhase;
				Bx[i][1] = b.x() * sinPhase;
				By[i][0] = b.y() * cosPhase;
				By[i][1] = b.y() * sinPhase;
				Bz[i][0] = b.z() * cosPhase;
				Bz[i][1] = b.z() * sinPhase;
			}
		}
	}

	// perform inverse FFT on each component
	fftw_execute(
			fftw_plan_dft_3d(n, n, n, Bx, Bx, FFTW_BACKWARD, FFTW_ESTIMATE));
	fftw_execute(
			fftw_plan_dft_3d(n, n, n, By, By, FFTW_BACKWARD, FFTW_ESTIMATE));
	fftw_execute(
			fftw_plan_dft_3d(n, n, n, Bz, Bz, FFTW_BACKWARD, FFTW_ESTIMATE));

	// normalize RMS of field to Brms
	double sumB2 = 0;
	for (unsigned int i = 0; i < n * n * n; i++)
		sumB2 += pow(Bx[i][0], 2) + pow(By[i][0], 2) + pow(Bz[i][0], 2);

	double weight = Brms / sqrt(sumB2 / (n * n * n));

	field.resize(n);
	for (int ix = 0; ix < n; ix++) {
		field[ix].resize(n);
		for (int iy = 0; iy < n; iy++) {
			field[ix][iy].resize(n);
		}
	}

	for (unsigned int ix = 0; ix < n; ix++)
		for (unsigned int iy = 0; iy < n; iy++)
			for (unsigned int iz = 0; iz < n; iz++) {
				int i = ix * n * n + iy * n + iz;
				field[ix][iy][iz] = Vector3(Bx[i][0], By[i][0], Bz[i][0])
						* weight;
			}

	fftw_free(Bx);
	fftw_free(By);
	fftw_free(Bz);
}

void periodicClamp(double pos, int n, int &lower, int &upper,
		double &frac2Lower, double &frac2Upper) {
	// closest lower and upper neighbour in a periodically continued grid
	lower = ((int(pos) % n) + n) % n;
	upper = (lower + 1) % n;
	frac2Lower = pos - floor(pos);
	frac2Upper = 1 - frac2Lower;
}

Vector3 TurbulentMagneticField::getField(const Vector3 &position) const {
	Vector3 r = (position - origin) / spacing;
	int ix, iX, iy, iY, iz, iZ;
	double fx, fX, fy, fY, fz, fZ;
	periodicClamp(r.x(), n, ix, iX, fx, fX);
	periodicClamp(r.y(), n, iy, iY, fy, fY);
	periodicClamp(r.z(), n, iz, iZ, fz, fZ);

	// trilinear interpolation
	// check: http://paulbourke.net/miscellaneous/interpolation/
	Vector3 b(0, 0, 0);
	//V000 (1 - x) (1 - y) (1 - z) +
	b += field[ix][iy][iz] * fX * fY * fZ;
	//V100 x (1 - y) (1 - z) +
	b += field[iX][iy][iz] * fx * fY * fZ;
	//V010 (1 - x) y (1 - z) +
	b += field[ix][iY][iz] * fX * fy * fZ;
	//V001 (1 - x) (1 - y) z +
	b += field[ix][iy][iZ] * fX * fY * fz;
	//V101 x (1 - y) z +
	b += field[iX][iy][iZ] * fx * fY * fz;
	//V011 (1 - x) y z +
	b += field[ix][iY][iZ] * fX * fy * fz;
	//V110 x y (1 - z) +
	b += field[iX][iY][iz] * fx * fy * fZ;
	//V111 x y z
	b += field[iX][iY][iZ] * fx * fy * fz;

	return b;
}

} // namespace mpc

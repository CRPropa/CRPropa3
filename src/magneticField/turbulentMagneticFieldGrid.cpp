#include "mpc/magneticField/turbulentMagneticFieldGrid.h"
#include "fftw3.h"

namespace mpc {

TurbulentMagneticFieldGrid::TurbulentMagneticFieldGrid(Vector3 origin,
		size_t samples, double spacing, double lMin, double lMax, double Brms,
		double powerSpectralIndex) :
		MagneticFieldGrid(origin, samples, spacing) {
	this->lMin = lMin;
	this->lMax = lMax;
	this->Brms = Brms;
	this->powerSpectralIndex = powerSpectralIndex;
	initialize();
}

void TurbulentMagneticFieldGrid::setSeed(int seed) {
	random.seed(seed);
	initialize();
}

void TurbulentMagneticFieldGrid::initialize() {
	size_t n = samples;

	// arrays to hold the complex vector components of the B-field
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
	double kMin = spacing / lMax;
	double kMax = spacing / lMin;
	Vector3 e1, e2, ek, b;
	Vector3 n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {

				i = ix * n * n + iy * n + iz;
				ek.set(K[ix], K[iy], K[iz]);
				k = ek.mag();

				if ((k < kMin) || (k > kMax)) {
					Bx[i][0] = 0;
					Bx[i][1] = 0;
					By[i][0] = 0;
					By[i][1] = 0;
					Bz[i][0] = 0;
					Bz[i][1] = 0;
					continue; // wave outside of turbulent range -> B(k) = 0
				}

				// construct an orthogonal base ek, e1, e2
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
				theta = 2 * M_PI * random.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// standard normal distributed amplitude weighted with k^alpha/2
				b *= random.randNorm() * pow(k, powerSpectralIndex / 2.);

				// uniform random phase
				phase = 2 * M_PI * random.rand();
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
	fftw_plan plan_x = fftw_plan_dft_3d(n, n, n, Bx, Bx, FFTW_BACKWARD,
			FFTW_ESTIMATE);
	fftw_execute(plan_x);
	fftw_destroy_plan(plan_x);

	fftw_plan plan_y = fftw_plan_dft_3d(n, n, n, By, By, FFTW_BACKWARD,
			FFTW_ESTIMATE);
	fftw_execute(plan_y);
	fftw_destroy_plan(plan_y);

	fftw_plan plan_z = fftw_plan_dft_3d(n, n, n, Bz, Bz, FFTW_BACKWARD,
			FFTW_ESTIMATE);
	fftw_execute(plan_z);
	fftw_destroy_plan(plan_z);

	// calculate normalization
	double sumB2 = 0;
	for (unsigned int i = 0; i < n * n * n; i++)
		sumB2 += pow(Bx[i][0], 2) + pow(By[i][0], 2) + pow(Bz[i][0], 2);
	double weight = Brms / sqrt(sumB2 / (n * n * n));

	// normalize and save real component to the grid
	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				int i = ix * n * n + iy * n + iz;
				grid[ix][iy][iz] = Vector3(Bx[i][0], By[i][0], Bz[i][0])
						* weight;
			}

	fftw_free(Bx);
	fftw_free(By);
	fftw_free(Bz);
}

double TurbulentMagneticFieldGrid::getRMSFieldStrength() const {
	return Brms;
}

double TurbulentMagneticFieldGrid::getPowerSpectralIndex() const {
	return powerSpectralIndex;
}

double TurbulentMagneticFieldGrid::getCorrelationLength() const {
	double r = lMin / lMax;
	double a = -powerSpectralIndex - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

} // namespace mpc

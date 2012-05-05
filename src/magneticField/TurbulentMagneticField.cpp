#include "mpc/magneticField/TurbulentMagneticField.h"
#include "fftw3.h"

namespace mpc {

void TurbulentMagneticField::setSeed(int seed) {
	random.seed(seed);
}

void TurbulentMagneticField::initialize(double lMin, double lMax, double Brms,
		double powerSpectralIndex) {
	this->lMin = lMin;
	this->lMax = lMax;
	this->Brms = Brms;
	this->powerSpectralIndex = powerSpectralIndex;

	size_t n = samples; // size of array
	size_t n2 = floor(n / 2) + 1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k, theta, phase, cosPhase, sinPhase;
	double kMin = spacing / lMax;
	double kMax = spacing / lMin;
	Vector3f b; // real b-field vector
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.set(K[ix], K[iy], K[iz]);
				k = ek.mag();

				// wave outside of turbulent range -> B(k) = 0
				if ((k < kMin) || (k > kMax)) {
					Bkx[i][0] = 0;
					Bkx[i][1] = 0;
					Bky[i][0] = 0;
					Bky[i][1] = 0;
					Bkz[i][0] = 0;
					Bkz[i][1] = 0;
					continue;
				}

				// construct an orthogonal base ek, e1, e2
				if (ek.isParallelTo(n0, 1e-6)) {
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

				Bkx[i][0] = b.x * cosPhase;
				Bkx[i][1] = b.x * sinPhase;
				Bky[i][0] = b.y * cosPhase;
				Bky[i][1] = b.y * sinPhase;
				Bkz[i][0] = b.z * cosPhase;
				Bkz[i][1] = b.z * sinPhase;
			}
		}
	}

	// in-place, complex to real, inverse Fourier transformation on each component
	// note that the last elements of B(x) are unused now
	float *Bx = (float*) Bkx;
	fftwf_plan plan_x = fftwf_plan_dft_c2r_3d(n, n, n, Bkx, Bx, FFTW_ESTIMATE);
	fftwf_execute(plan_x);
	fftwf_destroy_plan(plan_x);

	float *By = (float*) Bky;
	fftwf_plan plan_y = fftwf_plan_dft_c2r_3d(n, n, n, Bky, By, FFTW_ESTIMATE);
	fftwf_execute(plan_y);
	fftwf_destroy_plan(plan_y);

	float *Bz = (float*) Bkz;
	fftwf_plan plan_z = fftwf_plan_dft_c2r_3d(n, n, n, Bkz, Bz, FFTW_ESTIMATE);
	fftwf_execute(plan_z);
	fftwf_destroy_plan(plan_z);

	// calculate normalization
	double sumB2 = 0;
	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
				sumB2 += pow(Bx[i], 2) + pow(By[i], 2) + pow(Bz[i], 2);
			}
	double w = Brms / sqrt(sumB2 / (n * n * n));

	// normalize and save real component to the grid
	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
				Vector3f &b = get(ix, iy, iz);
				b.x = Bx[i] * w;
				b.y = By[i] * w;
				b.z = Bz[i] * w;
			}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);
}

double TurbulentMagneticField::getRMSFieldStrength() const {
	return Brms;
}

double TurbulentMagneticField::getPowerSpectralIndex() const {
	return powerSpectralIndex;
}

double TurbulentMagneticField::getCorrelationLength() const {
	double r = lMin / lMax;
	double a = -powerSpectralIndex - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

} // namespace mpc

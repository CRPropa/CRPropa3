#include "mpc/magneticField/TurbulentMagneticField.h"
#include "mpc/Units.h"

#include "fftw3.h"

namespace mpc {

TurbulentMagneticField::TurbulentMagneticField(Vector3d origin, double size,
		size_t samples) :
		MagneticFieldGrid(origin, size, samples) {
}

TurbulentMagneticField::TurbulentMagneticField(Vector3d origin, double size,
		size_t samples, double lMin, double lMax, double spectralIndex,
		double Brms) :
		MagneticFieldGrid(origin, size, samples) {
	this->lMin = lMin;
	this->lMax = lMax;
	this->spectralIndex = spectralIndex;
	initialize();
	normalize(Brms / getRMSFieldStrength());
}

void TurbulentMagneticField::setTurbulenceProperties(double lMin, double lMax,
		double spectralIndex) {
	this->lMin = lMin;
	this->lMax = lMax;
	this->spectralIndex = spectralIndex;
}

void TurbulentMagneticField::initialize(int seed) {
	random.seed(seed);
	initialize();
}

void TurbulentMagneticField::initialize() {
	if (lMin < 2 * spacing)
		throw std::runtime_error(
				"mpc::TurbulentMagneticField: lMin < 2 * spacing");
	if (lMin >= lMax)
		throw std::runtime_error("mpc::TurbulentMagneticField: lMin >= lMax");
	if (lMax > size / 2)
		throw std::runtime_error(
				"mpc::TurbulentMagneticField: lMax > size / 2");

	size_t n = samples; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

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
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getMag();

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
				if (ek.getAngleTo(n0) < 1e-3) {
					// ek parallel to (1,1,1)
					e1.setXYZ(-1., 1., 0);
					e2.setXYZ(1., 1., -2.);
				} else {
					// ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.getMag();
				e2 /= e2.getMag();

				// random orientation perpendicular to k
				theta = 2 * M_PI * random.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// standard normal distributed amplitude weighted with k^alpha/2
				b *= random.randNorm() * pow(k, spectralIndex / 2.);

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

	// save to grid
	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
				Vector3f &b = get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);
}

double TurbulentMagneticField::getPowerSpectralIndex() const {
	return spectralIndex;
}

double TurbulentMagneticField::getMinimumWavelength() const {
	return lMin;
}

double TurbulentMagneticField::getMaximumWavelength() const {
	return lMax;
}

double TurbulentMagneticField::getCorrelationLength() const {
	double r = lMin / lMax;
	double a = -spectralIndex - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

} // namespace mpc

#include "mpc/TurbulentMagneticField.h"
#include "mpc/MersenneTwister.h"
#include "fftw3.h"
#include <math.h>

namespace mpc {

TurbulentMagneticField::TurbulentMagneticField(Vector3 origin, size_t N,
		double spacing, double Brms, double spectralIndex, double Lmin,
		double Lmax) {
	this->origin = origin;
	this->N = N;
	this->spacing = spacing;
	this->Brms = Brms;
	this->spectralIndex = spectralIndex;
	this->kMin = 1. / Lmax;
	this->kMax = 1. / Lmin;
}

TurbulentMagneticField::~TurbulentMagneticField() {
}

void TurbulentMagneticField::initialize() {
	// random generator
	MTRand mtrand;

	// vector components of B-field
	fftw_complex *Bx, *By, *Bz;
	Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

	// N discrete possible wave numbers
	double K[N];
	for (int i = 0; i < N / 2; i++)
		K[i] = (float) i / N;
	for (int i = N / 2; i < N; i++)
		K[i] = (float) i / N - 1;

	int i;
	double k, theta, phase, cosPhase, sinPhase;
	Vector3 e1, e2, ek, b;
	Vector3 n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (int ix = 0; ix < N; ix++) {
		for (int iy = 0; iy < N; iy++) {
			for (int iz = 0; iz < N; iz++) {

				i = ix * N * N + iy * N + iz;
				ek.set(K[ix], K[iy], K[iz]);
				k = ek.mag();

				if ((k < kMin) || (k > kMax))
					continue; // wave outside of turbulent range -> B(k) = 0

				// construct an orthogonal base e1, e2, ek
				if ((ix == iy) && (iy == iz)) {
					// ek || (1,1,1)
					e1.set(-1., 1., 0);
					e2.set(1., 1., -2.);
				} else {
					// ek not || (1,1,1)
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
				fftw_complex a[2];
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
	fftw_plan p;
	p = fftw_plan_dft_3d(N, N, N, Bx, Bx, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	p = fftw_plan_dft_3d(N, N, N, By, By, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	p = fftw_plan_dft_3d(N, N, N, Bz, Bz, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	// normalize RMS of field to Brms
	double sumB2 = 0;
	for (unsigned int i = 0; i < N * N * N; i++)
		sumB2 += pow(Bx[i][0], 2) + pow(By[i][0], 2) + pow(Bz[i][0], 2);

	double weight = Brms / sqrt(sumB2 / (N * N * N));

	field.resize(N * N * N);
	for (unsigned int i = 0; i < N * N * N; i++) {
		field[i] = Vector3(Bx[i][0], By[i][0], Bz[i][0]) * weight;
	}
	fftw_free(Bx);
	fftw_free(By);
	fftw_free(Bz);
}

Vector3 TurbulentMagneticField::getField(const Vector3 &position) const {
	Vector3 r = (position - origin) / spacing;

	// index in a periodically continued grid
	int ix = ((int(r.x())%N)+N)%N;
	int iX = (ix + 1)%N;

	int iy = ((int(r.y())%N)+N)%N;
	int iY = (iy + 1)%N;

	int iz = ((int(r.z())%N)+N)%N;
	int iZ = (iz + 1)%N;

	double fx = r.x() - floor(r.x());
	double fX = 1 - fx;

	double fy = r.y() - floor(r.y());
	double fY = 1 - fy;

	double fz = r.z() - floor(r.y());
	double fZ = 1 - fz;

	Vector3 b(0.);
	int N2 = N * N;

	// check: http://paulbourke.net/miscellaneous/interpolation/
	//V000 (1 - x) (1 - y) (1 - z) +
	b += field[ix * N2 + iy * N + iz] * fX * fY * fZ;
	//V100 x (1 - y) (1 - z) +
	b += field[iX * N2 + iy * N + iz] * fx * fY * fZ;
	//V010 (1 - x) y (1 - z) +
	b += field[ix * N2 + iY * N + iz] * fX * fy * fZ;
	//V001 (1 - x) (1 - y) z +
	b += field[ix * N2 + iy * N + iZ] * fX * fY * fz;
	//V101 x (1 - y) z +
	b += field[iX * N2 + iy * N + iZ] * fx * fY * fz;
	//V011 (1 - x) y z +
	b += field[ix * N2 + iY * N + iZ] * fX * fy * fz;
	//V110 x y (1 - z) +
	b += field[iX * N2 + iY * N + iz] * fx * fy * fZ;
	//V111 x y z
	b += field[iX * N2 + iY * N + iZ] * fx * fy * fz;

	return b;
}

} // namespace mpc

#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#ifdef CRPROPA_HAVE_FFTW3F
#include "fftw3.h"

namespace crpropa {

SimpleGridTurbulence::SimpleGridTurbulence(const SimpleTurbulenceSpectrum &spectrum,
                                           const GridProperties &gridProp,
                                           unsigned int seed)
    : GridTurbulence(spectrum, gridProp, seed) {
	initTurbulence(gridPtr, spectrum.getBrms(), spectrum.getLmin(),
	               spectrum.getLmax(), -spectrum.getSindex() - 2, seed);
}

void SimpleGridTurbulence::initTurbulence(ref_ptr<Grid3f> grid, double Brms,
                                          double lMin, double lMax,
                                          double alpha, int seed) {
	
	Vector3d spacing = grid->getSpacing();

	checkGridRequirements(grid, lMin, lMax);

	size_t n = grid->getNx(); // size of array
	size_t n2 = (size_t)floor(n / 2) +
	            1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (size_t i = 0; i < n; i++)
		K[i] = (double)i / n - i / (n / 2);

	double kMin = spacing.x / lMax;
	double kMax = spacing.x / lMin;
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {
	
				Vector3f ek, e1, e2;  // orthogonal base

				size_t i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				double k = ek.getR();

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
				if (ek.isParallelTo(n0, float(1e-3))) {
					// ek parallel to (1,1,1)
					e1.setXYZ(-1., 1., 0);
					e2.setXYZ(1., 1., -2.);
				} else {
					// ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.getR();
				e2 /= e2.getR();

				// random orientation perpendicular to k
				double theta = 2 * M_PI * random.rand();
				Vector3f b = e1 * cos(theta) + e2 * sin(theta); // real b-field vector

				// normal distributed amplitude with mean = 0 and sigma =
				// k^alpha/2
				b *= random.randNorm() * pow(k, alpha / 2);

				// uniform random phase
				double phase = 2 * M_PI * random.rand();
				double cosPhase = cos(phase); // real part
				double sinPhase = sin(phase); // imaginary part

				Bkx[i][0] = b.x * cosPhase;
				Bkx[i][1] = b.x * sinPhase;
				Bky[i][0] = b.y * cosPhase;
				Bky[i][1] = b.y * sinPhase;
				Bkz[i][0] = b.z * cosPhase;
				Bkz[i][1] = b.z * sinPhase;
			} // for iz
		}     // for iy
	}         // for ix

	executeInverseFFTInplace(grid, Bkx, Bky, Bkz);

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}

} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

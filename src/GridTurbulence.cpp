#include "crpropa/GridTools.h"
#include "crpropa/GridTurbulence.h"
#include "crpropa/Random.h"
#include "crpropa/magneticField/MagneticField.h"

#include <map>

#ifdef CRPROPA_HAVE_FFTW3F
#include "fftw3.h"
#endif

namespace crpropa {

double turbulentCorrelationLength(double lMin, double lMax, double alpha) {
	double r = lMin / lMax;
	double a = -alpha - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

void checkGridRequirements(ref_ptr<VectorGrid> grid, double lMin, double lMax) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3d spacing = grid->getSpacing();

	if ((Nx != Ny) or (Ny != Nz))
		throw std::runtime_error("turbulentField: only cubic grid supported");
	if ((spacing.x != spacing.y) or (spacing.y != spacing.z))
		throw std::runtime_error("turbulentField: only equal spacing suported");
	if (lMin < 2 * spacing.x)
		throw std::runtime_error("turbulentField: lMin < 2 * spacing");
	if (lMin >= lMax)
		throw std::runtime_error("turbulentField: lMin >= lMax");
	if (lMax > Nx * spacing.x / 2)
		throw std::runtime_error("turbulentField: lMax > size / 2");
}

#ifdef CRPROPA_HAVE_FFTW3F

std::vector<std::pair<int, GridPrecision> > gridPowerSpectrum(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();

	double rms = rmsFieldStrength(grid);

	size_t n = Nx; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

	// construct the field in configuration space
	int i;
	double k;

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n);

	fftwf_complex *Bx = (fftwf_complex*) Bkx;
	fftwf_complex *By = (fftwf_complex*) Bky;
	fftwf_complex *Bz = (fftwf_complex*) Bkz;

	// save to temp
	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * n + iy * n + iz;
				Vector3<GridPrecision> &b = grid->get(ix, iy, iz);
				Bx[i][0] = b.x / rms;
				By[i][0] = b.y / rms;
				Bz[i][0] = b.z / rms;
			}
		}
	}

	// in-place, real to complex, inverse Fourier transformation on each component
	// note that the last elements of B(x) are unused now
	//fftwf_plan plan_x = fftwf_plan_dft_r2c_3d(n, n, n, Bx, Bkx, FFTW_ESTIMATE);
	fftwf_plan plan_x = fftwf_plan_dft_3d(n, n, n, Bx, Bkx, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(plan_x);
	fftwf_destroy_plan(plan_x);

	//fftwf_plan plan_y = fftwf_plan_dft_r2c_3d(n, n, n, By, Bky, FFTW_ESTIMATE);
	fftwf_plan plan_y = fftwf_plan_dft_3d(n, n, n, By, Bky, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(plan_y);
	fftwf_destroy_plan(plan_y);

	//fftwf_plan plan_z = fftwf_plan_dft_r2c_3d(n, n, n, Bz, Bkz, FFTW_ESTIMATE);
	fftwf_plan plan_z = fftwf_plan_dft_3d(n, n, n, Bz, Bkz, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(plan_z);
	fftwf_destroy_plan(plan_z);

	GridPrecision power;
	std::map<size_t, std::pair<GridPrecision, int> > spectrum;

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * n + iy * n + iz;
				k = static_cast<int>(std::floor(std::sqrt(ix*ix + iy*iy + iz*iz)));
				if (k > n/2. || k == 0)
					continue;
				power = ((Bkx[i][0]*Bkx[i][0] + Bkx[i][1]*Bkx[i][1]) +
					 (Bky[i][0]*Bky[i][0] + Bky[i][1]*Bky[i][1]) +
					 (Bkz[i][0]*Bkz[i][0] + Bkz[i][1]*Bkz[i][1]));
				if (spectrum.find(k) == spectrum.end()) {
					spectrum[k].first = power;
					spectrum[k].second = 1;
				}
				else {
					spectrum[k].first += power;
					spectrum[k].second += 1;
				}
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	std::vector<std::pair<int, GridPrecision> > points;
	for( std::map<size_t, std::pair<GridPrecision, int> >::iterator it = spectrum.begin(); it != spectrum.end(); ++it ) {
        	points.push_back(std::make_pair(
					it->first,
					(it->second).first/(it->second).second));
    	}

	return points;
}

void initHelicalTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, double alpha, int seed, double H) {

	checkGridRequirements(grid, lMin, lMax);

	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3d spacing = grid->getSpacing();

	size_t n = Nx; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k;

	// only used if there is a helicity
	double Bktot, Bkplus, Bkminus, thetaplus, thetaminus;

	double kMin = spacing.x / lMax;
	double kMax = spacing.x / lMin;
	Vector3f b; // real b-field vector
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getR();

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
				// (for helical fields together with the real transform the following convention
				// must be used: e1(-k) = e1(k), e2(-k) = - e2(k)
				if (ek.getAngleTo(n0) < 1e-3) { // ek parallel to (1,1,1)
					e1.setXYZ(-1, 1, 0);
					e2.setXYZ(1, 1, -2);
				} else { // ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.getR();
				e2 /= e2.getR();


				double Bkprefactor = mu0 / (4 * M_PI * pow(k, 3));
				Bktot = fabs(random.randNorm() * pow(k, alpha / 2));
				Bkplus  = Bkprefactor * sqrt((1 + H) / 2) * Bktot;
				Bkminus = Bkprefactor * sqrt((1 - H) / 2) * Bktot;
				thetaplus = 2 * M_PI * random.rand();
				thetaminus = 2 * M_PI * random.rand();
				double ctp = cos(thetaplus);
				double stp = sin(thetaplus);
				double ctm = cos(thetaminus);
				double stm = sin(thetaminus);

				Bkx[i][0] = ((Bkplus * ctp + Bkminus * ctm) * e1.x + (-Bkplus * stp + Bkminus * stm) * e2.x) / sqrt(2);
				Bkx[i][1] = ((Bkplus * stp + Bkminus * stm) * e1.x + ( Bkplus * ctp - Bkminus * ctm) * e2.x) / sqrt(2);
				Bky[i][0] = ((Bkplus * ctp + Bkminus * ctm) * e1.y + (-Bkplus * stp + Bkminus * stm) * e2.y) / sqrt(2);
				Bky[i][1] = ((Bkplus * stp + Bkminus * stm) * e1.y + ( Bkplus * ctp - Bkminus * ctm) * e2.y) / sqrt(2);
				Bkz[i][0] = ((Bkplus * ctp + Bkminus * ctm) * e1.z + (-Bkplus * stp + Bkminus * stm) * e2.z) / sqrt(2);
				Bkz[i][1] = ((Bkplus * stp + Bkminus * stm) * e1.z + ( Bkplus * ctp - Bkminus * ctm) * e2.z) / sqrt(2);

				Vector3f BkRe(Bkx[i][0], Bky[i][0], Bkz[i][0]);
				Vector3f BkIm(Bkx[i][1], Bky[i][1], Bkz[i][1]);
			} // for iz
		} // for iy
	} // for ix


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
				Vector3f &b = grid->get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}

void initTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, double alpha, int seed) {

	checkGridRequirements(grid, lMin, lMax);

	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3d spacing = grid->getSpacing();

	size_t n = Nx; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k;

	// parameters goes for non-helical calculations
	double theta, phase, cosPhase, sinPhase;

	double kMin = spacing.x / lMax;
	double kMax = spacing.x / lMin;
	Vector3f b; // real b-field vector
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getR();

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
				theta = 2 * M_PI * random.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// normal distributed amplitude with mean = 0 and sigma = k^alpha/2
				b *= random.randNorm() * pow(k, alpha / 2);

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
			} // for iz
		} // for iy
	} // for ix

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
				Vector3f &b = grid->get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}

void initTurbulenceWithBendover(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, double alpha, int seed, double lambda) {

	checkGridRequirements(grid, lMin, lMax);

	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3d spacing = grid->getSpacing();

	size_t n = Nx; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k;

	// parameters goes for non-helical calculations
	double theta, phase, cosPhase, sinPhase;

	double kMin = spacing.x / lMax;
	double kMax = spacing.x / lMin;
	Vector3f b; // real b-field vector
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getR();

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
				theta = 2 * M_PI * random.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// normal distributed amplitude with mean = 0 and sigma = k^alpha/2
				b *= random.randNorm() * k * lambda *
					pow(1+k*k*lambda*lambda, alpha / 4 - 1./2.);

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
			} // for iz
		} // for iy
	} // for ix

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
				Vector3f &b = grid->get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}
#endif // CRPROPA_HAVE_FFTW3F

} // namespace crpropa

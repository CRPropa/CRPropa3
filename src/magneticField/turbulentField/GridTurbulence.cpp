#include "crpropa/magneticField/turbulentField/GridTurbulence.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#ifdef CRPROPA_HAVE_FFTW3F

namespace crpropa {


GridTurbulence::GridTurbulence(const TurbulenceSpectrum &spectrum,
                               const GridProperties &gridProp,
                               unsigned int seed)
    : TurbulentField(spectrum), seed(seed) {
	initGrid(gridProp);
	checkGridRequirements(gridPtr, spectrum.getLmin(), spectrum.getLmax());
	initTurbulence();
}

void GridTurbulence::initGrid(const GridProperties &p) {
	gridPtr = new Grid3f(p);
}

Vector3d GridTurbulence::getField(const Vector3d &pos) const {
	return gridPtr->interpolate(pos);
}

const ref_ptr<Grid3f> &GridTurbulence::getGrid() const { return gridPtr; }

void GridTurbulence::initTurbulence() {

	Vector3d spacing = gridPtr->getSpacing();
	size_t n = gridPtr->getNx(); // size of array
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
		K[i] = ((double)i / n - i / (n / 2));

	// double kMin = 2*M_PI / lMax; // * 2 * spacing.x; // spacing.x / lMax;
	// double kMax = 2*M_PI / lMin; // * 2 * spacing.x; // spacing.x / lMin;
	double kMin = spacing.x / spectrum.getLmax();
	double kMax = spacing.x / spectrum.getLmin();
	auto lambda = 1 / spacing.x * 2 * M_PI;

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
				Vector3f b = e1 * std::cos(theta) + e2 * std::sin(theta); // real b-field vector

				// normal distributed amplitude with mean = 0
				b *= std::sqrt(spectrum.energySpectrum(k*lambda));
				
				// uniform random phase
				double phase = 2 * M_PI * random.rand();
				double cosPhase = std::cos(phase); // real part
				double sinPhase = std::sin(phase); // imaginary part

				Bkx[i][0] = b.x * cosPhase;
				Bkx[i][1] = b.x * sinPhase;
				Bky[i][0] = b.y * cosPhase;
				Bky[i][1] = b.y * sinPhase;
				Bkz[i][0] = b.z * cosPhase;
				Bkz[i][1] = b.z * sinPhase;
			} // for iz
		}     // for iy
	}         // for ix

	executeInverseFFTInplace(gridPtr, Bkx, Bky, Bkz);

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scaleGrid(gridPtr, spectrum.getBrms() /
	                       rmsFieldStrength(gridPtr)); // normalize to Brms
}

// Check the grid properties before the FFT procedure
void GridTurbulence::checkGridRequirements(ref_ptr<Grid3f> grid, double lMin,
                                           double lMax) {
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
	if (lMax > Nx * spacing.x) // before was (lMax > Nx * spacing.x / 2), why?
		throw std::runtime_error("turbulentField: lMax > size");
	if (lMax < lMin) 
		throw std::runtime_error("lMax < lMin");
}

// Execute inverse discrete FFT in-place for a 3D grid, from complex to real
// space
void GridTurbulence::executeInverseFFTInplace(ref_ptr<Grid3f> grid,
                                              fftwf_complex *Bkx,
                                              fftwf_complex *Bky,
                                              fftwf_complex *Bkz) {

	size_t n = grid->getNx(); // size of array
	size_t n2 = (size_t)floor(n / 2) +
	            1; // size array in z-direction in configuration space

	// in-place, complex to real, inverse Fourier transformation on each
	// component note that the last elements of B(x) are unused now
	float *Bx = (float *)Bkx;
	fftwf_plan plan_x = fftwf_plan_dft_c2r_3d(n, n, n, Bkx, Bx, FFTW_ESTIMATE);
	fftwf_execute(plan_x);
	fftwf_destroy_plan(plan_x);

	float *By = (float *)Bky;
	fftwf_plan plan_y = fftwf_plan_dft_c2r_3d(n, n, n, Bky, By, FFTW_ESTIMATE);
	fftwf_execute(plan_y);
	fftwf_destroy_plan(plan_y);

	float *Bz = (float *)Bkz;
	fftwf_plan plan_z = fftwf_plan_dft_c2r_3d(n, n, n, Bkz, Bz, FFTW_ESTIMATE);
	fftwf_execute(plan_z);
	fftwf_destroy_plan(plan_z);

	// save to grid
	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {
				size_t i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
				Vector3f &b = grid->get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}
}

Vector3f GridTurbulence::getMeanFieldVector() const {
	return meanFieldVector(gridPtr);
}

double GridTurbulence::getMeanFieldStrength() const {
	return meanFieldStrength(gridPtr);
}

double GridTurbulence::getRmsFieldStrength() const {
	return rmsFieldStrength(gridPtr);
}

std::array<float, 3> GridTurbulence::getRmsFieldStrengthPerAxis() const {
	return rmsFieldStrengthPerAxis(gridPtr);
}
	
std::vector<std::pair<int, float>> GridTurbulence::getPowerSpectrum() const {
	return gridPowerSpectrum(gridPtr);
}

void GridTurbulence::dumpToFile(std::string filename) const {
	dumpGrid(gridPtr, filename);
}

} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

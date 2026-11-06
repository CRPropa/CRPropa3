#include "crpropa/magneticField/UF26Field.h"
#include "crpropa/GridTools.h"
#include "crpropa/Units.h"
#ifdef CRPROPA_HAVE_FFTW3F
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"
#endif

#include "kiss/logger.h"

#include <cmath>
#include <stdexcept>

// local helper functions and constants
namespace uf26 {
	// Internal units used in this code. Conversion to CRPropa units
	// happens only in the public get...() methods.
	const double kpc = 1;
	const double microgauss = 1;
	// galactocentric radius of the Sun (Gravity Collaboration 2019)
	const double kSunRadius = 8.178 * kpc;
	// arcosh(e): sech(arcosh(e) z/z_d) = 1/e at z = z_d
	const double kArcoshE = std::acosh(std::exp(1.));

	// Gaussian G(x; sigma) of Eq. (16)
	inline double gauss(const double x, const double sigma) {
		return std::exp(-0.5 * x * x / (sigma * sigma));
	}
}

namespace crpropa {

UF26Field::UF26Field(const ModelType mt) : fModelType(mt) {
	// common to both models
	fDiskRmax = 18 * uf26::kpc;
	fDiskWmax = 2 * uf26::kpc;

	switch (fModelType) {
	case expDisk:
		fDiskB = 4.4 * uf26::microgauss;
		fDiskZ = 1.0 * uf26::kpc;
		fDiskL = 15 * uf26::kpc;
		fDiskRc = 1 * uf26::kpc;
		fHasAnnulus = false;
		break;
	case ringDisk:
		fDiskB = 4.1 * uf26::microgauss;
		fDiskZ = 0.9 * uf26::kpc;
		fDiskL = 100 * uf26::kpc;
		fDiskRc = 0;
		fAnnulusB = 4.1 * uf26::microgauss;
		fAnnulusR = 4.5 * uf26::kpc;
		fAnnulusW = 2.0 * uf26::kpc;
		fAnnulusZ = 1.9 * uf26::kpc;
		fHasAnnulus = true;
		break;
	default:
		throw std::runtime_error("UF26Field: unknown model type");
	}

	fSunRadiusSmoothed =
		std::sqrt(uf26::kSunRadius * uf26::kSunRadius + fDiskRc * fDiskRc);
}

#ifdef CRPROPA_HAVE_FFTW3F
void UF26Field::randomTurbulent(int seed) {
	// Kolmogorov spectrum with B_rms = 1 (scaled by b(x)) between 4 and
	// 136 pc, i.e. Lc = 30 pc, on a periodic grid of 256^3 points with
	// 2 pc spacing
	const double lMin = 4 * parsec;
	const double lMax = 136 * parsec;
	const double Brms = 1;
	const double spacing = 2 * parsec;
	const size_t gridN = 256;

	auto spectrum = SimpleTurbulenceSpectrum(Brms, lMin, lMax);
	auto gp = GridProperties(Vector3d(0.), gridN, spacing);
	auto tf = SimpleGridTurbulence(spectrum, gp, seed);
	setTurbulentGrid(tf.getGrid());
}
#endif

void UF26Field::setTurbulentGrid(ref_ptr<Grid3f> grid) {
	if (!grid.valid())
		throw std::runtime_error("UF26Field::setTurbulentGrid: invalid grid");
	const double rms = rmsFieldStrength(grid);
	if (std::fabs(rms - 1) > 0.01) {
		KISS_LOG_WARNING << "UF26Field: turbulent grid has Brms = " << rms
			<< ", expected 1.";
	}
	turbulentGrid = grid;
}

ref_ptr<Grid3f> UF26Field::getTurbulentGrid() {
	return turbulentGrid;
}

bool UF26Field::isUsingTurbulentField() const {
	return turbulentGrid.valid();
}

UF26Field::ModelType UF26Field::getModelType() const {
	return fModelType;
}

std::string UF26Field::getModelName() const {
	return fModelType == expDisk ? "expDisk" : "ringDisk";
}

Vector3d UF26Field::getField(const Vector3d &pos) const {
	return getTurbulentField(pos);
}

Vector3d UF26Field::getTurbulentField(const Vector3d &pos) const {
	if (!turbulentGrid.valid())
		throw std::runtime_error("UF26Field: no turbulent grid set, "
			"run e.g. randomTurbulent()");
	const double b = getTurbulentStrength(pos);
	// skip the grid interpolation where the model vanishes
	if (b == 0)
		return Vector3d(0.);
	return Vector3d(turbulentGrid->interpolate(pos)) * b;
}

double UF26Field::getTurbulentStrength(const Vector3d &pos) const {
	const Vector3d p = pos / kpc * uf26::kpc;
	// rms strengths of independent components add in quadrature
	return std::hypot(rmsDisk(p), rmsAnnulus(p)) / uf26::microgauss * microgauss;
}

double UF26Field::rmsDisk(const Vector3d &pos) const {
	const double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
	// E_d(r)/E_d(R_sun), Eq. (13)
	const double radial =
		std::exp(-(std::sqrt(r * r + fDiskRc * fDiskRc) - fSunRadiusSmoothed) / fDiskL);
	// outer cutoff S_d(r)
	const double cutoff = 1 / (1 + std::exp((r - fDiskRmax) / fDiskWmax));
	// vertical profile V_d(z), Eq. (14)
	const double vertical = 1 / std::cosh(uf26::kArcoshE * pos.z / fDiskZ);
	double b = fDiskB * radial * cutoff * vertical;
	// ringDisk: disk is phased out toward the center following the inner
	// flank of the annulus, G(min(r - r_a, 0); w_a)
	if (fHasAnnulus && r < fAnnulusR)
		b *= uf26::gauss(r - fAnnulusR, fAnnulusW);
	return b;
}

double UF26Field::rmsAnnulus(const Vector3d &pos) const {
	if (!fHasAnnulus)
		return 0;
	const double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
	// vertical Gaussian with 1/e height z_a, i.e. G(z; z_a/sqrt(2))
	const double zz = pos.z / fAnnulusZ;
	return fAnnulusB * uf26::gauss(r - fAnnulusR, fAnnulusW) * std::exp(-zz * zz);
}

} // namespace crpropa

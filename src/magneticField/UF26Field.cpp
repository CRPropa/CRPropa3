#include "crpropa/magneticField/UF26Field.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/Units.h"

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

UF26Field::UF26Field(const ModelType mt, ref_ptr<MagneticField> turbulence) :
	UF26Field(mt) {
	setTurbulentField(turbulence);
}

void UF26Field::setTurbulentField(ref_ptr<MagneticField> turbulence) {
	const TurbulentField* tf = dynamic_cast<const TurbulentField*>(turbulence.get());
	if (!tf)
		throw std::runtime_error("UF26Field::setTurbulentField: field is not a "
			"TurbulentField, use setTurbulentField(field, brms)");
	// read Brms now: TurbulentField only holds a reference to its spectrum
	setTurbulentField(turbulence, tf->getBrms());
}

void UF26Field::setTurbulentField(ref_ptr<MagneticField> field, double brms) {
	if (!field.valid())
		throw std::runtime_error("UF26Field::setTurbulentField: invalid field");
	if (!(brms > 0))
		throw std::runtime_error("UF26Field::setTurbulentField: brms must be > 0");
	fTurbulence = field;
	fTurbulenceBrms = brms;
}

ref_ptr<MagneticField> UF26Field::getTurbulentField() const {
	return fTurbulence;
}

bool UF26Field::hasTurbulentField() const {
	return fTurbulence.valid();
}

UF26Field::ModelType UF26Field::getModelType() const {
	return fModelType;
}

std::string UF26Field::getModelName() const {
	return fModelType == expDisk ? "expDisk" : "ringDisk";
}

Vector3d UF26Field::getField(const Vector3d &pos) const {
	if (!fTurbulence.valid())
		throw std::runtime_error("UF26Field::getField: no turbulent field set, "
			"see setTurbulentField()");
	const double b = getRmsField(pos);
	// skip the (possibly expensive) turbulence where the model vanishes
	if (b == 0)
		return Vector3d(0.);
	return fTurbulence->getField(pos) * (b / fTurbulenceBrms);
}

double UF26Field::getRmsField(const Vector3d &pos) const {
	const Vector3d p = pos / kpc * uf26::kpc;
	// rms strengths of independent components add in quadrature
	return std::hypot(rmsDisk(p), rmsAnnulus(p)) / uf26::microgauss * microgauss;
}

double UF26Field::getDiskRmsField(const Vector3d &pos) const {
	return rmsDisk(pos / kpc * uf26::kpc) / uf26::microgauss * microgauss;
}

double UF26Field::getAnnularRmsField(const Vector3d &pos) const {
	return rmsAnnulus(pos / kpc * uf26::kpc) / uf26::microgauss * microgauss;
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

#include "crpropa/magneticField/TF17Field.h"
#include "crpropa/Units.h"

#include <algorithm>

namespace crpropa {

TF17Field::TF17Field() {
	// disk parameters
	useDiskField = true;
	a_disk = 0.031 / kpc / kpc;
	r1_disk = 3 * kpc;
	B1_disk = 32.0 * muG;
	phi_star_disk = -31 * M_PI / 180;
	H_disk = 0.054 * kpc;

	// halo parameters
	useHaloField = true;
	a_halo = 0.33 / kpc / kpc;
	z1_halo = 0 * kpc;
	B1_halo = 9.0 * muG;
	L_halo = 1.2 * kpc;
	phi_star_halo = 198 * M_PI / 180;
	m = 1; // for bisymmetric halo model

	// shared parameters
	p_0 = -9.1 * M_PI / 180;
	cot_p0 = cos(p_0) / sin(p_0);
	H_p = 1.2 * kpc;
	L_p = 50 * kpc;		// best fit simply says L_p > 40 kpc:  L_p ≫ L_e can hardly be constrained by the FDobs map.
}

void TF17Field::setUseDiskField(bool use) {
	useDiskField = use;
}

void TF17Field::setUseHaloField(bool use) {
	useHaloField = use;
}

bool TF17Field::isUsingDiskField() {
	return useDiskField;
}

bool TF17Field::isUsingHaloField() {
	return useHaloField;
}

double TF17Field::zscale(const double& z) const {
	return 1 + (fabs(z) / H_p) * (fabs(z) / H_p);
}

double TF17Field::shiftedWindingFunction(const double& r, const double& z) const {
	return cot_p0 * log(1 - exp(-r / L_p)) / zscale(z);
}

double TF17Field::radialFieldScale(const double& B1, const double& r1, const double& z1, const double& phi1) const {
	// This term occures is parameterizations of models A and B
	return B1 * exp(-fabs(z1) / H_disk) * cos(m * (phi1 - shiftedWindingFunction(r1, z1) - phi_star_disk));
}

double TF17Field::verticalFieldScale(const double& B1, const double& r1, const double& z1, const double& phi1) const {
	// This term occures is parameterizations of models C and D
	return B1 * exp(-r1 / L_halo) * cos(m * (phi1 - shiftedWindingFunction(r1, z1) - phi_star_halo));
}

double TF17Field::azimuthalFieldComponent(const double& r, const double& z, const double& B_r, const double& B_z) const {
	double r_ = r / L_p;
	double B_phi = cot_p0 / zscale(z) * r_ * exp(-r_) / (1 - exp(-r_)) * B_r - \
								 2 * cot_p0 * z / (H_p * H_p) / (zscale(z) * zscale(z)) * r * log(1 - exp(-r_)) * B_z;
	return B_phi;
}

Vector3d TF17Field::getHaloField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const {
	// Model C1, bisymmetric (m=1)
	Vector3d b(0.);

	double r1_halo = r / (1 + a_halo * z * z);
	double phi1_halo = phi - shiftedWindingFunction(r, z) + shiftedWindingFunction(r1_halo, z1_halo);
	// B components in (r, phi, z)
	double r_ = r > std::numeric_limits<double>::epsilon() ? (r1_halo / r) : 1;		// avoid zero division
	double B_r = 2 * a_halo * r_ * r_ * r1_halo * z * verticalFieldScale(B1_halo, r1_halo, z1_halo, phi1_halo);
	double B_z = r_ * r_ * verticalFieldScale(B1_halo, r1_halo, z1_halo, phi1_halo);
	double B_phi = azimuthalFieldComponent(r, z, B_r, B_z);
	// Convert to (x, y, z) components
	b.x = - (B_r * cosPhi - B_phi * sinPhi);	// flip x-component at the end
	b.y = B_r * sinPhi + B_phi * cosPhi;
	b.z = B_z;
	return b;
}

Vector3d TF17Field::getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const {
	// Model Ad1
	Vector3d b(0.);
	double B_r = 0;
	double B_phi = 0;
	double B_z = 0;

	double z1_disk = z * (1 + a_disk * r1_disk * r1_disk) / (1 + a_disk * r * r);
	if (r > r1_disk) {
		double phi1_disk = phi - shiftedWindingFunction(r, z) + shiftedWindingFunction(r1_disk, z1_disk);
		// B components in (r, phi, z)
		double z_ = fabs(z) > std::numeric_limits<double>::epsilon() ? (z1_disk / z) : 1;		// avoid zero division
		double B_r0 = radialFieldScale(B1_disk, r1_disk, z1_disk, phi1_disk);
		B_r = (r1_disk / r) * z_ * B_r0;
		B_z = 2 * a_disk * r1_disk * z1_disk / (1+ a_disk * r * r) * B_r0;
		B_phi = azimuthalFieldComponent(r, z, B_r, B_z);
	} else {
		// within r = 3 kpc, the field lines are straight in direction g_phi + phi_star_disk
		double g_phi = shiftedWindingFunction(r1_disk, z1_disk);
		double B_amp = B1_disk * exp(-fabs(z1_disk) / H_disk);
		B_r = cos(g_phi + phi_star_disk) * B_amp;
		B_phi = sin(g_phi + phi_star_disk) * B_amp;
	}
	// Convert to (x, y, z) components
	b.x = - (B_r * cosPhi - B_phi * sinPhi);	// flip x-component at the end
	b.y = B_r * sinPhi + B_phi * cosPhi;
	b.z = B_z;
	return b;
}

Vector3d TF17Field::getField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
	double phi = M_PI - pos.getPhi(); // azimuth in our convention
	// double cosPhi = pos.x / r;
	double cosPhi = cos(phi);
	// double sinPhi = pos.y / r;
	double sinPhi = sin(phi);

	Vector3d b(0.);
	if (useDiskField)
		b += getDiskField(r, pos.z, phi, sinPhi, cosPhi);	// disk field
	if (useHaloField)
		b += getHaloField(r, pos.z, phi, sinPhi, cosPhi);	// halo field
	return b;
}

} // namespace crpropa

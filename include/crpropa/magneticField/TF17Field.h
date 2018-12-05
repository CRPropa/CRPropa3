#ifndef CRPROPA_TF17FIELD_H
#define CRPROPA_TF17FIELD_H

#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

/**
 Implements the TF2017 galactic magnetic field model, consisting of large-scale regular disk and halo fields.

 Currently only best fit values of the field parameters are implemented: Halo model C1; disk model Ad1.
 The field is defined in the usual galactocentric coordinate system with the Galactic center at the origin,
 the x-axis pointing in the opposite direction of  the Sun, and the z-axis pointing towards Galactic north.

 See: Terral, Ferriere 2017 - Constraints from Faraday rotation on the magnetic field structure in the galactic halo, DOI: 10.1051/0004-6361/201629572, arXiv:1611.10222
 */

class TF17Field: public MagneticField {
private:

	// disk parameters
	bool useDiskAd1;
	double a_disk;
	double r1_disk;
	double B1_disk;
	double phi_star_halo;
	double H_disk;

	// halo parameters
	bool useHaloC;
	double a_halo;
	double z1_halo;
	double cot_p0;
	double B1_halo;
	double L_halo;
	double phi_star_halo;
	double m;

	// universal parameters
	double p_0;
	double H_p;
	double L_p;

	void SetParams();

public:
	TF17Field();

	double zscale(const double& z) const;
	double shiftedWindingFunction(const double& r, const double& z) const;
	double radialFieldScale(const double& B1, const double& r1, const double& z1, const double& phi1) const;
	double verticalFieldScale(const double& B1, const double& r1, const double& z1, const double& phi1) const;
	double azimuthalFieldComponent(const double& r, const double& z, const double& B_r, const double& B_z) const;

	Vector3d getHaloField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;
	Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;
	Vector3d getField(const Vector3d& pos) const;
};

} // namespace crpropa

#endif // CRPROPA_TF17FIELD_H

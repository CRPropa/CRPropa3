#ifndef CRPROPA_TF17FIELD_H
#define CRPROPA_TF17FIELD_H

#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
using namespace std;

/**
 Implements the TF2017 galactic magnetic field model, consisting of large-scale regular disk and halo fields.
 The field is defined in the usual galactocentric coordinate system with the Galactic center at the origin,
 the x-axis pointing in the opposite direction of  the Sun, and the z-axis pointing towards Galactic north.
 See: Terral, Ferriere 2017 - Constraints from Faraday rotation on the magnetic field structure in the galactic halo, DOI: 10.1051/0004-6361/201629572, arXiv:1611.10222

 */

class TF17Field: public MagneticField {
private:
    // model definition
    bool C0;
    bool C1;
    bool Ad1;
    bool Bd1;
    bool Dd1;

	// disk parameters
	bool useDiskField;
	double a_disk;
	double z1_disk;
	double r1_disk;
	double B1_disk;
	double L_disk;
	double phi_star_disk;
	double H_disk;

	// halo parameters
	bool useHaloField;
	double a_halo;
	double z1_halo;
	double cot_p0;
	double B1_halo;
	double L_halo;
	double phi_star_halo;

	// universal parameters
	double p_0;
	double H_p;
	double L_p;

	void SetParams();

public:
	TF17Field(string halo_model="C1", string disk_model="Ad1");

	void setUseDiskField(bool use);
	void setUseHaloField(bool use);
	bool isUsingDiskField();
	bool isUsingHaloField();

	Vector3d getField(const Vector3d& pos) const;
	Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;
	Vector3d getHaloField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;

	double azimuthalFieldComponent(const double& r, const double& z, const double& B_r, const double& B_z) const;
	double radialFieldScale(const double& B1, const double& r1, const double& z1, const double& phi1) const;
	double verticalFieldScale(const double& B1, const double& r1, const double& z1, const double& phi1, const double& L, const int& m) const;
	double shiftedWindingFunction(const double& r, const double& z) const;
	double zscale(const double& z) const;
};

} // CRPROPA NAMESPACE

#endif // CRPROPA_TF17FIELD_H

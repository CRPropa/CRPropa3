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

enum class TF17DiskModel {Ad1, Bd1, Dd1};
enum class TF17HaloModel {C0, C1};


class TF17Field: public MagneticField {
private:
    TF17DiskModel disk_model; 
    TF17HaloModel halo_model;

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
	double B1_halo;
	double L_halo;
	double phi_star_halo;

	// universal parameters
	double p_0;
	double cot_p0;
	double H_p;
	double L_p;

	// security to avoid 0 division
	double epsilon;

	void SetParams();

public:
	TF17Field(TF17DiskModel disk_model_=TF17DiskModel::Ad1, TF17HaloModel halo_model_=TF17HaloModel::C1);

	void setUseDiskField(bool use);
	void setUseHaloField(bool use);
	bool isUsingDiskField();
	bool isUsingHaloField();

    void set_B1_disk(const double B1);
    void set_z1_disk(const double z1);
    void set_r1_disk(const double r1);
    void set_H_disk(const double H);
    void set_L_disk(const double L);
    void set_a_disk(const double a);
    void set_phi_star_disk(const double phi);

    void set_B1_halo(const double B1);
    void set_z1_halo(const double z1);
    void set_L_halo(const double L);
    void set_a_halo(const double a);
    void set_phi_star_halo(const double phi);

    void set_Hp(const double H);
    void set_Lp(const double L);
    void set_p0(const double p0);

    string getDiskModel() const;
    string getHaloModel() const;

	Vector3d getField(const Vector3d& pos) const;
	Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;
	Vector3d getHaloField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;

	double azimuthalFieldComponent(const double& r, const double& z, const double& B_r, const double& B_z) const;
	double radialFieldScale(const double& B1, const double& phi_star, const double& z1, const double& phi, const double& r, const double& z) const;
	double shiftedWindingFunction(const double& r, const double& z) const;
	double zscale(const double& z) const;
};

} // CRPROPA NAMESPACE

#endif // CRPROPA_TF17FIELD_H

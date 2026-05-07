#ifndef CRPROPA_TF17FIELD_H
#define CRPROPA_TF17FIELD_H

#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
using namespace std;
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class TF17Field
 @brief TF17Field galactic magnetic field model
 
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
    /** Constructor
        @param disk_model_     model to use for the disk (model are defined in enum class TF17DiskModel: Ad1, Bd1, Dd1)
        @param halo_model_     model to use for the halo (model are defined in enum class TF17HaloModel: C0, C1)
     */
	TF17Field(TF17DiskModel disk_model_=TF17DiskModel::Ad1, TF17HaloModel halo_model_=TF17HaloModel::C1);

	void setUseDiskField(bool use);
	void setUseHaloField(bool use);
	bool isUsingDiskField();
	bool isUsingHaloField();

    /**@brief   set reference amplitude for the disk magnetic field
     *          Automatically set at initialization
     * @param B1    magnetic field amplitude in Gauss
     */
    void set_B1_disk(const double B1);
    /**@brief   set reference height for the disk magnetic field (model Dd1) 
     *          Automatically set at initialization
     * @param z1    height in kpc 
     */
    void set_z1_disk(const double z1);
    /**@brief   set reference radius for the disk magnetic field (models Ad1 and Bd1)
     *          Automatically set at initialization
     * @param r1    radius in kpc 
     */
    void set_r1_disk(const double r1);
    /**@brief   set scale height of Br for the disk magnetic field (models Ad1 and Bd1)
     *          Automatically set at initialization
     * @param H    height in kpc 
     */
    void set_H_disk(const double H);
    /**@brief   set scale length of Bz for the disk magnetic field (models Dd1)
     *          Automatically set at initialization
     * @param L    length in kpc 
     */
    void set_L_disk(const double L);
    /**@brief   set opening parameter for poloidal lines for the disk magnetic field (models Ad1)
     *          Automatically set at initialization
     * @param a     opening parameter in / kpc / kpc
     */
    void set_a_disk(const double a);
    /**@brief   set orientation angle of the azimuthal pattern for the disk magnetic field
     *          Automatically set at initialization
     * @param phi     opening parameter in rad
     */
    void set_phi_star_disk(const double phi);

    /**@brief   set reference amplitude for the halo magnetic field
     *          Automatically set at initialization
     * @param B1    magnetic field amplitude in Gauss
     */
    void set_B1_halo(const double B1);
    /**@brief   set reference height for the halo magnetic field
     *          Automatically set at initialization
     * @param z1    height in kpc 
     */
    void set_z1_halo(const double z1);
    /**@brief   set scale length of Bz for the halo magnetic field
     *          Automatically set at initialization
     * @param L    length in kpc 
     */
    void set_L_halo(const double L);
    /**@brief   set opening parameter for poloidal lines for the halo magnetic field
     *          Automatically set at initialization
     * @param a     opening parameter in / kpc / kpc
     */
    void set_a_halo(const double a);
    /**@brief   set orientation angle of the azimuthal pattern for the disk magnetic field
     *          Automatically set at initialization
     * @param phi     opening parameter in rad
     */
    void set_phi_star_halo(const double phi);

    /**@brief   set pitch angle origin (all models)
     *          Automatically set at initialization
     * @param p0    height in rad
     */
    void set_p0(const double p0);
    /**@brief   set scale height of winding function (all models)
     *          Automatically set at initialization
     * @param H    height in kpc 
     */
    void set_Hp(const double H);
    /**@brief   set scale length of winding function (all models)
     *          Automatically set at initialization
     * @param L    length in kpc 
     */
    void set_Lp(const double L);

    /**@brief   Get the disk model used as a string
     * @return  disk model
     */
    string getDiskModel() const;
    /**@brief   Get the halo model used as a string
     * @return  halo model
     */
    string getHaloModel() const;

	Vector3d getField(const Vector3d& pos) const;
	Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;
	Vector3d getHaloField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;

    /**@brief   Compute the azimuthal field component Bphi as define by equation 28 in TF17
     * @param r     radius in cylindrical coordinates
     * @param z     radius in cylindrical coordinates
     * @param B_r    radial component of the magnetic field at position (r,z)
     * @param B_z    height component of the magnetic field at position (r,z)
     * @return  the value of the azimuthal field component Bphi
     */
	double azimuthalFieldComponent(const double& r, const double& z, const double& B_r, const double& B_z) const;

    /**@brief   Compute the scaling of the disk magnetic field amplitude (equation 30, models Ad1
     * and Bd1)
     * @param B1        magnetic field amplitude in Gauss
     * @param phi_star  orientation angle of the azimuthal pattern in rad
     * @param z1        scale height in kpc
     * @param phi       cylindrical coordinates
     * @param r         cylindrical coordinates
     * @param z         cylindrical coordinates
     * @return  the radial field component Br(r1, z1, phi1)
     */
	double radialFieldScale(const double& B1, const double& phi_star, const double& z1, const double& phi, const double& r, const double& z) const;

    /**@brief   Compute the shifted winding function (equation 23)
     * @param r         cylindrical coordinates
     * @param z         cylindrical coordinates
     * @return  g_phi, the shifted winding function
     */
	double shiftedWindingFunction(const double& r, const double& z) const;

    /**@brief   Compute the height scaling appearing numerous times in the equations
     * @param z         cylindrical coordinates
     * @return  z_scale = 1 + (z/Hp)^2
     */
	double zscale(const double& z) const;
};
/**@}*/
} // CRPROPA NAMESPACE

#endif // CRPROPA_TF17FIELD_H

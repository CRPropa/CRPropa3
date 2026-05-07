#ifndef CRPROPA_JF12FIELD_H
#define CRPROPA_JF12FIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"
#include "kiss/logger.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class JF12Field
 @brief JF12Field galactic magnetic field model

 Implements the JF2012 magnetic field model, consisting of a large-scale regular
 and random (striated) field and a small-scale random (turbulent) field.
 See:
 Jansson 2012a, ApJ. 757, A New Model of the Galactic Magnetic Field
 Jansson 2012b, arXiv:1210.7820, The Galactic Magnetic Field

 All three components may individually turned on and off.
 Currently only best fit values of the field paramaters are implemented and
 cannot be changed.

 The field is defined in the usual galactocentric coordinate system with the
 Galactic center at the origin, the x-axis pointing in the opposite direction of
 the Sun, and the z-axis pointing towards Galactic north.

 The regular field components (disk, toroidal halo and polodial halo field)
 may be turned on and off individually.
 */
class JF12Field: public MagneticField {
protected:
	bool useRegularField;
	bool useStriatedField;
	bool useTurbulentField;
	bool useDiskField;
	bool useToroidalHaloField;
	bool useXField;

	// disk spiral arms
	double rArms[8];       // radii where each arm crosses the negative x-axis
	double pitch;          // pitch angle
	double sinPitch, cosPitch, tanPitch, cotPitch, tan90MinusPitch;

	// Regular field ----------------------------------------------------------
	// disk
	double bDisk[11];      // field strengths of the 8 arms at r=5 kpc; additional entries added for periodic closure in JF12FieldSolenoidal
	double bRing;          // ring field strength 3<r<5 kpc
	double hDisk, wDisk;   // disk/halo transistion and width
	// toroidal halo
	double bNorth, bSouth; // northern, southern halo field strength
	double rNorth, rSouth; // northern, southern transistion radius
	double wHalo, z0;      // transistion width and vertical scale height
	// poloidal halo
	double bX;             // field strength at origin
	double thetaX0;        // constant elevation angle at r > rXc, z = 0
	double sinThetaX0, cosThetaX0, tanThetaX0, cotThetaX0;
	double rXc;            // radius of varying elevation angle region
	double rX;             // exponential scale height

	// Striated field ---------------------------------------------------------
	double sqrtbeta;       // relative strength of striated field
	ref_ptr<Grid1f> striatedGrid;

	// Turbulent field --------------------------------------------------------
	ref_ptr<Grid3f> turbulentGrid;
	// disk
	double bDiskTurb[8]; // field strengths in arms at r=5 kpc
	double bDiskTurb5;   // field strength at r<5kpc
	double zDiskTurb;	 // Gaussian scale height of disk
	// halo
	double bHaloTurb; // halo field strength
	double rHaloTurb; // exponential scale length
	double zHaloTurb; // Gaussian scale height

public:
	JF12Field();

	// Create and set a random realization for the striated field
	void randomStriated(int seed = 0);

#ifdef CRPROPA_HAVE_FFTW3F
	// Create a random realization for the turbulent field
	void randomTurbulent(int seed = 0);
#endif

	/**
	 * Set a striated grid and activate the striated field component
	 * @param grid	scalar grid containing random +1/-1 values, 100 parsec grid spacing
	 */
	void setStriatedGrid(ref_ptr<Grid1f> grid);

	/**
	 * Set a turbulent grid and activate the turbulent field component
	 * @param grid	vector grid containing a random field of Brms = 1
	 */
	void setTurbulentGrid(ref_ptr<Grid3f> grid);

	ref_ptr<Grid1f> getStriatedGrid();
	ref_ptr<Grid3f> getTurbulentGrid();

	void setUseRegularField(bool use);
	virtual void setUseStriatedField(bool use);
	virtual void setUseTurbulentField(bool use);
	void setUseDiskField(bool use);
	void setUseToroidalHaloField(bool use);
	void setUseXField(bool use);

	bool isUsingRegularField();
	bool isUsingStriatedField();
	bool isUsingTurbulentField();
	bool isUsingDiskField();
	bool isUsingToroidalHaloField();
	bool isUsingXField();

	double logisticFunction(const double& x, const double& x0, const double& w) const;

	// Regular field components
	Vector3d getRegularField(const Vector3d& pos) const;
	virtual Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;
	Vector3d getToroidalHaloField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const;
	virtual Vector3d getXField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const;

	// Regular and striated field component
	Vector3d getStriatedField(const Vector3d& pos) const;

	// Brms of the turbulent field
	double getTurbulentStrength(const Vector3d& pos) const;

	// Turbulent field component
	Vector3d getTurbulentField(const Vector3d& pos) const;

	// All set field components
	Vector3d getField(const Vector3d& pos) const;
};
/** @} */


/**
 * \addtogroup MagneticFields
 * @{
 */
/**
 @class PlanckJF12bField
 @brief PlanckJF12bField: the JF12 galactic magnetic field model with corrections
 suggested by the Planck Collaboration

 See: Planck Collaboration, "Planck intermediate results. XLII.
 Large-scale Galactic magnetic fields", A&A 596 (2016) A103;
 arXiv:1601.00546

 This variant of the JF12 field uses only different parameters compared to the standard JF12 implementation.
  */
class PlanckJF12bField: public JF12Field {
	public:
		PlanckJF12bField();
};
/** @} */

} // namespace crpropa

#endif // CRPROPA_JF12FIELD_H

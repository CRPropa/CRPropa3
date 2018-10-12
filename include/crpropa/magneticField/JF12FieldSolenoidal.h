#ifndef CRPROPA_JF12FIELDSOLENOIDAL_H
#define CRPROPA_JF12FIELDSOLENOIDAL_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"
#include <crpropa/Units.h>

namespace crpropa {

/**
 @class JF12FieldSolenoidal
 @brief JF12FieldSolenoidal galactic magnetic field model
 
 Implements a modified JF2012 magnetic field model, see documentation of JF12Field for basics.
 
 A solenoidal ring-spiral transition (of width delta = 3 kpc per default) for the disk field
 and a smooth X-field at z = 0 (parabolic field lines for abs(z) < zs; zs = 500 parsec per default)
 were added in order to smooth the field and avoid violations of magnetic flux conservation, see:
 Kleimann et al, 2018, arXiv:1809.07528, Solenoidal improvements for the JF12 Galactic magnetic field model.
 
 The regular field components (disk, toroidal halo and polodial halo field) 
 may be turned on and off individually for tests.
 
 The turbulent field component is the exact same as in the initial JF12 implementation
 and should be used with care since the new regular and old turbulent field
 do not match in the ring-spiral transition region.
 */
 
class JF12FieldSolenoidal: public MagneticField {
private:

	// allow for enabling and disabling of the individual components
	bool useRegular;
	bool useStriated;
	bool useTurbulent;
	
	bool useDisk;
	bool useToroidalHalo;
	bool useX;
	
	// height paramter of the modified X-field
	double zS;
	
	// disk spiral arms
	double rArms[8]; // radii where each arm crosses the negative x-axis
	
	// angles at which the dividing spirals intersect the r1-ring, 
	// longer array due to cyclicity
	double phi0Arms[11];
	
	// lower bound for phi integration
	double phi0;
	// corresponding index in phi0Arms
	int idx0;
	
	// phi integration H(phi) = phiCoeff[j] + bDisk[j] * phi
	double phiCoeff[10]; // for phi0Arms[j-1] < phi < phi0Arms[j], hence only 10
	double corr; // correction term enforcing H(phi0) = 0
	
	// inner and outer boundaries of disk field
	double r1; 
	double r2;
	
	// transitions region boundaries of disk field
	double r1s;
	double r2s;
	
	double pitch;          // pitch angle
	double sinPitch, cosPitch, tanPitch, cotPitch, tan90MinusPitch;

	// Regular field ----------------------------------------------------------
	// disk
	double bDisk[11];       // field strengths of arms at r=5 kpc
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
	ref_ptr<ScalarGrid> striatedGrid;

	// Turbulent field --------------------------------------------------------
	ref_ptr<VectorGrid> turbulentGrid;
	// disk
	double bDiskTurb[8]; // field strengths in arms at r=5 kpc
	double bDiskTurb5;   // field strength at r<5kpc
	double zDiskTurb;	 // Gaussian scale height of disk
	// halo
	double bHaloTurb; // halo field strength
	double rHaloTurb; // exponential scale length
	double zHaloTurb; // Gaussian scale height

public:

	JF12FieldSolenoidal(double delta=(3*kpc), double zs=(0.5*kpc));

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
	void setStriatedGrid(ref_ptr<ScalarGrid> grid);

	/**
	 * Set a turbulent grid and activate the turbulent field component
	 * @param grid	vector grid containing a random field of Brms = 1
	 */
	void setTurbulentGrid(ref_ptr<VectorGrid> grid);

	ref_ptr<ScalarGrid> getStriatedGrid();
	ref_ptr<VectorGrid> getTurbulentGrid();

	void setUseRegular(bool use);
	void setUseStriated(bool use);
	void setUseTurbulent(bool use);
	void setUseDisk(bool use);
	void setUseToroidalHalo(bool use);
	void setUseX(bool use);
	
	void setDelta(double d);
	void setZs(double z);

	bool isUsingRegular();
	bool isUsingStriated();
	bool isUsingTurbulent();
	bool isUsingDisk();
	bool isUsingToroidalHalo();
	bool isUsingX();
	
	double getDelta() const;
	double getZs() const;
	
	// continue spiral field to r=20kpc without transition at the outer boundary
	// you may reactivate the transition afterwards with setDelta()
	void deactivateOuterTransition();
	
	// scaling function for disk and toroidal halo field; same as in initial JF12
	double logisticFunction(const double x, const double x0, const double w) const;

	// Regular field component
	Vector3d getRegularField(const Vector3d& pos) const;

	// Regular and striated field component
	Vector3d getStriatedField(const Vector3d& pos) const;

	// Brms of the turbulent field
	double getTurbulentStrength(const Vector3d& pos) const;

	// Turbulent field component
	Vector3d getTurbulentField(const Vector3d& pos) const;

	// Total field
	Vector3d getField(const Vector3d& pos) const;
	
	// Individual regular field components
	Vector3d getXField(const double r, const double z, const double sinPhi, const double cosPhi) const;
	Vector3d getDiskField(const double r, const double z, const double phi, const double sinPhi, const double cosPhi) const;
	Vector3d getToroidalHaloField(const double r, const double z, const double sinPhi, const double cosPhi) const;
	
	// transition polynomial p_delta(r)
	double p(const double r) const;
	
	// transition polynomial derivative
	double q(const double r) const;
	
	// evaluate H-Integral
	double PhiIntegralH(const double r, const double phi) const;
	
	// return correct field strength b_j for current spiral arm
	double getSpiralStrength(const double r, const double phi) const;
	
};

} // namespace crpropa

#endif // CRPROPA_JF12FIELDSOLENOIDAL_H

#ifndef MPC_JF2012FIELD_H
#define MPC_JF2012FIELD_H

#include "mpc/magneticField/MagneticField.h"
#include "mpc/Grid.h"

namespace mpc {

class JF2012Field: public MagneticField {
	bool useStriated;
	bool useTurbulent;

	// disk spiral arms
	double rArms[8];       // radii where each arm crosses the negative x-axis
	double pitch;          // pitch angle
	double sinPitch, cosPitch, tan90MinusPitch;

	// Regular field ---------------------------------------------------------
	// disk
	double bDisk[8];       // field strengths of arms at r=5 kpc
	double bRing;          // ring field strength 3<r<5 kpc
	double hDisk, wDisk;   // disk/halo transistion and width
	// toroidal halo
	double bNorth, bSouth; // northern, southern halo field strength
	double rNorth, rSouth; // northern, southern transistion radius
	double wHalo, z0;      // transistion width and vertical scale height
	// poloidal halo
	double bX;             // field strength at origin
	double thetaX0;        // constant elevation angle at r > rXc, z = 0
	double sinThetaX0, cosThetaX0, tanThetaX0;
	double rXc;            // radius of varying elevation angle region
	double rX;             // exponential scale height

	// Striated field --------------------------------------------------------
	double sqrtbeta;       // relative strength of striated field
	ref_ptr<ScalarGrid> striatedGrid;

	// Small-scale turbulent field -------------------------------------------
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
	JF2012Field();

	void randomStriated(int seed = 0);
	void randomTurbulent(int seed = 0);

	void setStriatedGrid(ref_ptr<ScalarGrid> grid);
	void setTurbulentGrid(ref_ptr<VectorGrid> grid);

	ref_ptr<ScalarGrid> getStriatedGrid();
	ref_ptr<VectorGrid> getTurbulentGrid();

	Vector3d getRegularField(const Vector3d& pos) const;
	Vector3d getStriatedField(const Vector3d& pos) const;
	Vector3d getTurbulentField(const Vector3d& pos) const;

	Vector3d getField(const Vector3d& pos) const;
};

} // namespace mpc

#endif // MPC_JF2012FIELD_H

#ifndef MPC_JF2012FIELD_H
#define MPC_JF2012FIELD_H

#include "mpc/magneticField/MagneticField.h"
#include "mpc/Grid.h"

namespace mpc {

class JF2012Field: public MagneticField {
	// disk field
	double bDisk[8]; // field strengths of arms at r=5 kpc, b8 is determined from other 7
	double rArms[8];     // radii where each arm crosses the negative x-axis
	double pitch;        // pitch angle
	double sinPitch, cosPitch, tan90MinusPitch;
	double bRing;        // ring field strength 3<r<5 kpc
	double hDisk, wDisk; // disk/halo transistion and width

	// toroidal halo
	double bNorth, bSouth; // northern, southern halo field strength
	double rNorth, rSouth; // northern, southern transistion radius
	double wHalo, z0;      // transistion width and vertical scale height

	// X halo
	double bX;           // field strength at origin
	double thetaX0;      // constant elevation angle at r > rXc, z = 0
	double sinThetaX0, cosThetaX0, tanThetaX0;
	double rXc;          // radius of varying elevation angle region
	double rX;           // exponential scale height

public:
	JF2012Field();
	Vector3d getField(const Vector3d& pos) const;
};

class JF2012TurbulentField: public MagneticField {
	double b5;
	double bDisk[8];	// field strengths of arms at r=5 kpc
	double rArms[8];	// radii where each arm crosses the negative x-axis
	double pitch;
	double tan90MinusPitch;
	double r0;	// radial attenuation length
	double z0;	// vertical attenuation length
	double z0S; // vertical attenuation length for spiral field
	ref_ptr<VectorGrid> grid; // turbulent field grid

public:
	JF2012TurbulentField();
	double getModulation(const Vector3d& pos) const;
	Vector3d getField(const Vector3d& pos) const;
};

}

#endif // MPC_JF2012FIELD_H

#ifndef CRPROPA_PSHIRKOVFIELD_H
#define CRPROPA_PSHIRKOVFIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {

/**
 @class PshirkovField
 @brief PshirkovField galactic magnetic field model

 Implements the Pshirkov magnetic field model, consisting of a large-scale regular
 (disk and halo) field.
 See:
 Pshirkov 2011, arXiv:1103.0814, Deriving global structure of the Galactic Magnetic Field from Faraday Rotation Measures of extragalactic sources

 Two models, the Axisymmetric (ASS) and the Bisymmetric (BSS) disk model can be chosen.
 Currently only best fit values of the field paramaters are implemented and
 cannot be changed.

 The field is defined in the usual galactocentric coordinate system with the
 Galactic center at the origin, the x-axis pointing in the opposite direction of
 the Sun, and the z-axis pointing towards Galactic north.
 */

class PshirkovField: public MagneticField {
private:
	bool useASS;
	bool useBSS;
	bool useBHMHalo;

	// disk parameters
	double pitch_ASS;	// pitch angle
	double pitch_BSS;
	double b_ASS, b_BSS, d;	// b = 1/tan(p) // d means the distance to the first field reversal
	double R_sun;		// distance of the Sun to the galactic center
	double R_c;		// radius of central region
	double theta_ASS, theta_BSS;		// theta = b*ln(1+d/R_sun) - pi/2
	double z0;		// vertical thickness of the MF in the galactic disk
	double B0;		// scale of the magnetic field strength
	
	// halo (north)
	double z0_n;		// vertical position of the halo field
	double R0_n;		// radial position of the halo field
	double B0_n;		// scale of the magnetic field strength in the halo
	double z11_n;		// vertical thickness towards the galactic disc
	double z12_n;		// vertical thickness off the galactic disk

	// halo (south)
	double z0_s;		// SAME with south halo
	double R0_s;
	double B0_s_ASS;
	double B0_s_BSS;
	double z11_s;
	double z12_s;
	
public:
	PshirkovField();

	void setUseASS(bool use);
	void setUseBSS(bool use);
	void setUseBHMHalo(bool use);

	bool isUsingASS();
	bool isUsingBSS();
	bool isUsingBHMHalo();
	
	double getSouthernFieldstrength() const;

	// Regular ASS field
	Vector3d getASSField(const Vector3d& pos) const;

	// Regular BSS field
	Vector3d getBSSField(const Vector3d& pos) const;
	
	// Regular BHM Halo field (Sun et al. 2008)
	Vector3d getBHMHaloField(const Vector3d& pos) const;
	
	// All set field components
	Vector3d getField(const Vector3d& pos) const;
};

} // namespace crpropa

#endif // CRPROPA_PSHIRKOVFIELD_H

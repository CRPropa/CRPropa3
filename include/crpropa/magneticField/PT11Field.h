#ifndef CRPROPA_PSHIRKOVFIELD_H
#define CRPROPA_PSHIRKOVFIELD_H

#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class PT11Field
 @brief Implements Pshirkov2011 galactic magnetic field model

 Implements the Pshirkov2011 galactic magnetic field model, consisting of large-scale regular disk and halo fields.
 For the disk field an axisymmetric (ASS) and the bisymmetric (BSS, default) model can be chosen.
 For the halo field the BHM Halo field model (Sun et al. 2008) is considered.

 Currently only best fit values of the field parameters are implemented.
 The field is defined in the usual galactocentric coordinate system with the Galactic center at the origin, the x-axis pointing in the opposite direction of  the Sun, and the z-axis pointing towards Galactic north.

 See: Pshirkov, Tinyakov, Kronberg Newton-McGee 2011 - Deriving global structure of the Galactic Magnetic Field from Faraday Rotation Measures of extragalactic sources, DOI: 10.1088/0004-637X/738/2/192, arXiv:1103.0814
 */

class PT11Field: public MagneticField {
private:
	bool useASS;  // switch for axisymmetric spiral field (ASS)
	bool useBSS;  // switch for bisymmetric spiral field (BSS)
	bool useHalo; // switch for halo field

	// disk parameters
	double pitch, cos_pitch, sin_pitch, PHI, cos_PHI;  // pitch angle parameters
	double d;     // distance to first field reversal
	double R_sun; // distance between sun and galactic center
	double R_c;   // radius of central region
	double z0_D;    // vertical thickness in the galactic disk
	double B0_D;    // magnetic field scale

	// halo parameters
	double z0_H;  // halo vertical position
	double R0_H;  // halo radial position
	double B0_Hn; // halo magnetic field scale (north)
	double B0_Hs; // halo magnetic field scale (south)
	double z11_H; // halo vertical thickness towards disc
	double z12_H; // halo vertical thickness off the disk

	void SetParams();

public:
	PT11Field();

	void setUseASS(bool use);
	void setUseBSS(bool use);
	void setUseHalo(bool use);

	bool isUsingASS();
	bool isUsingBSS();
	bool isUsingHalo();

	Vector3d getField(const Vector3d& pos) const;
};
/**@}*/
} // namespace crpropa

#endif // CRPROPA_PSHIRKOVFIELD_H

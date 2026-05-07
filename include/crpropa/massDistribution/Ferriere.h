#ifndef CRPROPA_FERRIERE_H
#define CRPROPA_FERRIERE_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <string>

namespace crpropa {
/**
 @class Ferriere
 @brief model of the distribution of hydrogen in the Milky Way
  Here in model Ferriere 2007
  seperated in 2 regions (inner, outer). The border is for R=3 kpc in galactocentric radius.
  model is discribed in
outer: ApJ, 497, 759
inner:	arxiv:	astro-ph/0702532
*/
class Ferriere: public Density {
private:
	// standard for all types of distribution
	bool isforHI = true;
	bool isforHII = true;
	bool isforH2 = true;
	double Rsun = 8.5 * kpc;  // distance sun-galactic center

public:
	/** Coordinate transformation for the CentralMolecularZone region. Rotation arround z-axis such that X is the major axis and Y is the minor axis
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return position in local coordinates for the CMZ region
	*/
	Vector3d CMZTransformation(const Vector3d &position) const;
	
	/** Coordinate transformation for the galactic bulge disk region in galactic center. Rotation arround the x-axis, the y'-axis and the x''-axis. Difened with X along the major axis, Y along the minor axis and Z along the northern normal
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return position in local coordinates for the GB disk region
	*/
	Vector3d DiskTransformation(const Vector3d &position) const;

	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density in parts/m^3, only acitvated parts are summed up */
	double getDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of atomic hydrogen in parts/m^3 */
	double getHIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of ionised hydrogen in parts/m^3 */
	double getHIIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of molecular hydrogen in parts/m^3 */
	double getH2Density(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return nucleon density in parts/m^3, only activated parts are summed up and H2 is weighted twice */
	double getNucleonDensity(const Vector3d &position) const;

	/** changes activation status for atomic hydrogen */
	void setIsForHI(bool HI);
	/** changes activation status for ionised hydrogen */
	void setIsForHII(bool HII);
	/** changes activation status for molecular hydrogen */
	void setIsForH2(bool H2);

	/** @return activation status for atomic hydrogen */
	bool getIsForHI();
	/** @return activation status for ionised hydrogen */
	bool getIsForHII();
	/** @return activation status for molecular hydrogen */
	bool getIsForH2();

	std::string getDescription();
};

}  // namespace crpropa

#endif  // CRPROPA_FERRIERE_H



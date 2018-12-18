#ifndef CRPROPA_NAKANISHI_H
#define CRPROPA_NAKANISHI_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <string>

namespace crpropa {
/**
 @class Nakanishi
 @brief zylindrical symetrical model of the density distribution of the Milkyway for atomic (HI) and molecular (H2) hydrogen
 	Modell for HI arXiv:astro-ph/0304338
	Modell for H2 arxiv:astro-ph/0610769
	fit of the models given in arXiv:1607.07886
*/
class Nakanishi: public Density{
private:
	bool isforHI = true;
	bool isforHII = false;
	bool isforH2 = true;

public:
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density in parts/m^3, only acitvated parts are summed up */
	double getDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of atomic hydrogen in parts/m^3 */
	double getHIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of molecular hydrogen in parts/m^3 */
	double getH2Density(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return nucleon density in parts/m^3, only activated parts are summed up and H2 is weighted twice */
	double getNucleonDensity(const Vector3d &position) const;

	/** the scaleheight over the galactic plane of atomic hydrogen is fitted by polynome of degree 3
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return scaleheight at given position */
	double getHIScaleheight(const Vector3d &position)const;
	/** the plane density is fittet by two exponential components with e^-R and e^-(R^2)
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return plane density in parts/m^3 */
	double getHIPlanedensity(const Vector3d &position)const;

	/** the scaleheight over the galactic plane of molecular hydrogen is fitted by exponential function
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return scaleheight at given position */
	double getH2Scaleheight(const Vector3d &position)const;
	/** the plane density is fittet by two exponential components
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return plane density in parts/m^3 */
	double getH2Planedensity(const Vector3d &position)const;

	/** changes activation status for atomic hydrogen */
	void setIsForHI(bool HI);
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

#endif  // CRPROPA_NAKANISHI_H




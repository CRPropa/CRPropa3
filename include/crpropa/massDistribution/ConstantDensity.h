#ifndef CRPROPA_CONSTANTDENSITY_H
#define CRPROPA_CONSTANTDENSITY_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <sstream>
#include <string>

namespace crpropa {

/*
@class ConstantDensity
@brief Density module for Constant densities in HI, HII and H2 component.
*/
class ConstantDensity: public Density {
private:
	// default mode: all density types set to 0 and no activ component
	double HIdensitynumber  = 0;  /**< density for atomic hydrogen */
	double HIIdensitynumber = 0;  /**< density for ioniesd hydrogen */
	double H2densitynumber  = 0;  /**< density for molecular hydrogen */

	bool isHI = false;  /**< If true, HI is used for sum up in getDensity */
	bool isHII = false;  /**< If true, HII is used for sum up in getDensity */
	bool isH2 = false;  /**< If true, H2 is used for sum up in getDensity */

public:
	/** Constructor for constant density
	 @param HI density for atomic hydrogen
	 @param HII density for ionised hydrogen
	 @param H2 density for molecular hydrogen
	 */
	ConstantDensity(double HI, double HII, double H2);

	/** Get density at a given position.
	 @param position 	position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density in parts/m^3, sum up all activated parts
	 */
	double getDensity(const Vector3d &position) const;
	/** Get HI density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns (constant) density of HI in parts/m^3
	 */
	double getHIDensity(const Vector3d &position) const;
	/** Get HII density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns (constant) density of HII in parts/m^3 
	 */
	double getHIIDensity(const Vector3d &position) const;
	/** Get H2 density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns (constant) density of H2 in parts/m^3 
	*/
	double getH2Density(const Vector3d &position) const;
	/** Get density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns number of nucleons/m^3, sum up all activated parts and weights H2 twice 
	 */
	double getNucleonDensity(const Vector3d &position) const;

	/** Status of HI -- active or not.
	 @returns Boolean flag with activation status of HI 
	 */
	bool getIsForHI();
	/** Status of HII -- active or not.
	 @returns Boolean flag with activation status of HII
	 */
	bool getIsForHII();
	/** Status of H2 -- active or not.
	 @returns Boolean flag with activation status of H2
	 */
	bool getIsForH2();

	/** Change HI status and the value of the density.
	 @param activate 		new activation status
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHI(bool activate, double densityNumber);
	/** Change HI status and keep density unaltered.
	 @param activate 		new activation status
	 */
	void setHI(bool activate);
	/** Change HI density and keep activation status unaltered
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHI(double densityNumber);

	/** Change HII status and the value of the density.
	 @param activate 		new activation status
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHII(bool activate, double densityNumber);
	/** Change HII status and keep density unaltered.
	 @param activate 		new activation status
	 */
	void setHII(bool activate);
	/** Change HII density and keep activation status unaltered
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHII(double densityNumber);

	/** Change H2 status and the value of the density.
	 @param activate 		new activation status
	 @param densityNumber	new density [in units of 1/meter ^ 3]
 	 */
	void setH2(bool activate, double densityNumber);
	/** Change H2 status and keep density unaltered.
	 @param activate 		new activation status
	 */
	void setH2(bool activate);
	/** Change H2 density and keep activation status unaltered
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setH2(double densityNumber);

	std::string getDescription();
};

}  // namespace crpropa

#endif  // CRPROPA_CONSTANTDENSITY_H



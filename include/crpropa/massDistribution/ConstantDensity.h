#ifndef CRPROPA_CONSTANTDENSITY_H
#define CRPROPA_CONSTANTDENSITY_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <sstream>

#include "kiss/logger.h"

namespace crpropa {

/*
@class ConstantDensity 
@brief Density module for Constant densitys in HI, HII and H2 component. 
*/
class ConstantDensity: public Density {

private:
	// default mode: all density types set to 0 and no activ component
	double HIdensitynumber  = 0;	/**< density for atomic hydrogen */
	double HIIdensitynumber = 0;	/**< density for ioniesd hydrogen */
	double H2densitynumber  = 0;	/**< density for molecular hydrogen */
	
	bool isHI = false;	/**< If true, HI is used for sum up in getDensity */
	bool isHII = false;	/**< If true, HII is used for sum up in getDensity */
	bool isH2 = false;	/**< If true, H2 is used for sum up in getDensity */

public:
	/** Constructor for constant density
	 @param HI density for atomic hydrogen
	 @param HII density for ionised hydrogen 
	 @param H2 density for molecular hydrogen
	 */
	ConstantDensity(double HI, double HII, double H2);
	
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0) 
	 @return density in parts/m^3, sum up all activated parts */
	double getDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0) 
	 @return (constant) density of HI in parts/m^3 */
	double getHIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0) 
	 @return (constant) density of HII in parts/m^3 */
	double getHIIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0) 
	 @return (constant) density of H2 in parts/m^3 */
	double getH2Density(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0) 
	@return number of nucleons/m^3, sum up all activated parts and weights H2 twice */
	double getNucleonDensity(const Vector3d &position) const;
	
	/** @return activation status of HI */
	bool getIsForHI();
	/** @return activation status of HII */
	bool getIsForHII();
	/** @return activation status of H2 */
	bool getIsForH2();
	
	/** change HI status and density number 
	 @param activate new activation status
	 @param densitynumber new densitynumber	*/
	void setHI(bool activate, double densitynumber);
	/** change HI status and keep density number as it is 
	 @param activate new activation status	*/
	void setHI(bool activate);
	/** change HI density number and keep status as it is 
	@param densitynumber new densitynumber */
	void setHI(double densitynumber);
	
	/** change HII status and density number 
	@param activate new activation status */
	void setHII(bool activate, double densitynumber);
	/** change HII status and keep dansity number as it is 
	@param activate new activation status */
	void setHII(bool activate);
	/** change HII dendisty number and keep status as it is 
	@param densitynumber new densitynumber */
	void setHII(double densitynumber);
	
	/** change H2 status and denstiy number 
	@param activate new activation status */
	void setH2(bool activate, double densitynumber);
	/** change H2 status and keep density number as it is 
	@param activate new activation status */
	void setH2(bool activate);
	/** change H2 density number and keep status as it is
	@param densitynumber new densitynumber */
	void setH2(double densitynumber);
	
	std::string getDescription();
};

} //namespace

#endif //CRPROPA_CONSTANTDENSITY_H



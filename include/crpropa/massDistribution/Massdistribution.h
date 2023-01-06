#ifndef CRPROPA_MASSDISTRIBUTION_H
#define CRPROPA_MASSDISTRIBUTION_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Grid.h"

#include "kiss/logger.h"

#include <vector>

namespace crpropa {

/**
 @class DensityList
 @brief Superposition of density models.
 The addDensity function adds a new density to the list.
 The getDensity function handles the activated types in loaded densities, whereas get(type)Density disregards the activation state.
*/
class DensityList: public Density {
private:
	std::vector<ref_ptr<Density> > DensityList ;

public:
	/** Add new density to list.
	 @param density density to add
	*/
	void addDensity(ref_ptr<Density> density);

	/** Get density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density in parts/m^3, sum up densities from added densities 
	*/
	double getDensity(const Vector3d &position) const;
	/** Get HI density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density of HI at given position in parts/m^3, sum up all HI densities from added densities
	 */
	double getHIDensity(const Vector3d &position) const;
	/** Get HII density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density of HII at given position in parts/m^3, sum up all HII densities from added densities */
	double getHIIDensity(const Vector3d &position) const;
	/** Get H2 density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density of H2 at given position in parts/m^3, sum up all H2 densities from added densities 
	 */
	double getH2Density(const Vector3d &position) const;
	/** Get the density of nucleons.
	 This is the number of nucleons per volume, summed up all activated density and weight molecular hydrogyen twice
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density of nucleons at given position in parts/m^3, sum up all nucleon densities from added densities 
	 */
	double getNucleonDensity(const Vector3d &position) const;
};

class DensityGrid: public Density {
private: 
	ref_ptr<Grid1f> grid; //< Grid with data
	bool isForHI, isForHII, isForH2; 
	void checkAndWarn(); //< raise a warning if all density types are deactivated.

public:
	DensityGrid(ref_ptr<Grid1f> grid, bool isForHI = false, bool isForHII = false, bool isForH2 = false);

	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;
	double getDensity(const Vector3d &position) const;
	double getNucleonDensity(const Vector3d &position) const;

	bool getIsForHI();
	bool getIsForHII();
	bool getIsForH2();

	void setIsForHI(bool b);
	void setIsForHII(bool b);
	void setIsForH2(bool b);
	
	void setGrid(ref_ptr<Grid1f> grid);
};

}  // namespace crpropa

#endif  // CRPROPA_MASSDISTRIBUTION_H
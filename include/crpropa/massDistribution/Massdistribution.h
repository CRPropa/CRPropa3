#ifndef CRPROPA_MASSDISTRIBUTION_H
#define CRPROPA_MASSDISTRIBUTION_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/Vector3.h"

#include "kiss/logger.h"

#include <vector>

namespace crpropa {
/**
 @class DensityList
 @brief Superposition of density models.
 the addDensity function adds a new density to the list.
 The getDensity function cares about acitvated types in loaded densitys. The get(typ)Density doesn't care.
*/
class DensityList: public Density {
private:
	std::vector<ref_ptr<Density> > DensityList ;

public:
	/** add new density to list
	@param density density to add*/
	void addDensity(ref_ptr<Density> density);

	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density in parts/m^3, sum up densities from added densities */
	double getDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of HI at given position in parts/m^3, sum up all HI densities from added densities */
	double getHIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of HII at given position in parts/m^3, sum up all HII densities from added densities */
	double getHIIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of H2 at given position in parts/m^3, sum up all H2 densities from added densities */
	double getH2Density(const Vector3d &position) const;
	/** NucleonDensity is the number of nucleons per Volume, sum up all activated density and weight molecular hydrogyen twice
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of nucleons at given position in parts/m^3, sum up all nucleon densities from added densities */
	double getNucleonDensity(const Vector3d &position) const;
};

}  // namespace crpropa

#endif  // CRPROPA_MASSDISTRIBUTION_H



#ifndef CRPROPA_CORDES_H
#define CRPROPA_CORDES_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <string>

namespace crpropa {
/**
	@class Cordes
	@brief Cylindrical symetrical model of the density of ionised hydrogen (HII) of the Milky Way
	Cordes et al., 1991, Nature 353,737
	*/
class Cordes: public Density {
private:
	// DO NOT CHANGE model type!
	bool isforHI = false;
	bool isforHII = true;
	bool isforH2 = false;

public:
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density in parts/m^3 */
	double getDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of ionised hydrogen in parts/m^3, equal getDensity thus no other type is included for Cordes */
	double getHIIDensity(const Vector3d &position) const;
	/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 @return density of nucleons in parts/m^3, equal getDensity thus only HII is included for Cordes */
	double getNucleonDensity(const Vector3d &position) const;

	/** @return activation status of HI */
	bool getIsForHI();
	/** @return activation status of HII */
	bool getIsForHII();
	/** @return activation status of H2 */
	bool getIsForH2();

	std::string getDescription();
};

}  // namespace crpropa

#endif  // CRPROPA_CORDES_H

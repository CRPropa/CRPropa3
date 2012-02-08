#ifndef TURBULENTMAGNETICFIELD_H_
#define TURBULENTMAGNETICFIELD_H_

#include "mpc/magneticField/magneticFieldGrid.h"

namespace mpc {

/**
 @class TurbulentMagneticFieldGrid
 @brief Random turbulent magnetic field on a cubic grid with trilinear interpolation.

 This class creates a random magnetic field with the following properties. \n
 mean and rms strength: \n
 \f$ <\vec{B}> = 0 \f$ \n
 \f$ <B^2> = B_{rms}^2\f$ \n
 turbulent range and spectrum: \n
 \f$ \vec{B}(\vec{k}) \sim k^{-11/3} \f$ for \f$ k_{min} < k < k_{max} \f$ \n
 \f$ \vec{B}(\vec{k}) = 0 \f$ otherwise \n
 solenoidity: \n
 \f$ \nabla \vec{B}(\vec{x}) = 0 \f$ \n
 isotropy: \n
 same coherence length in all directions
 */
class TurbulentMagneticFieldGrid: public MagneticFieldGrid {
public:
	TurbulentMagneticFieldGrid(Vector3 origin, size_t samples, double spacing,
			double Brms, double lMin, double lMax, double powerSpectralIndex,
			int seed);
protected:
	void initialize();
	double Brms, lMin, lMax, powerSpectralIndex;
	int seed;
};

} // namespace mpc

#endif /* TURBULENTMAGNETICFIELD_H_ */

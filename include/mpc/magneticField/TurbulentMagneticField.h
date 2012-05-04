#ifndef MPC_TURBULENTMAGNETICFIELD_H_
#define MPC_TURBULENTMAGNETICFIELD_H_

#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/Random.h"

namespace mpc {

/**
 @class TurbulentMagneticField
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
class TurbulentMagneticField: public MagneticFieldGrid {
public:
	TurbulentMagneticField(Vector3d origin, size_t samples, double spacing) :
			MagneticFieldGrid(origin, samples, spacing) {
	}
	void initialize(double lMin, double lMax, double Brms,
			double powerSpectralIndex = -11. / 3.);
	void setSeed(int seed);
	double getRMSFieldStrength() const;
	double getCorrelationLength() const;
	double getPowerSpectralIndex() const;

protected:
	double Brms, lMin, lMax, powerSpectralIndex;
	Random random;
};

} // namespace mpc

#endif /* MPC_TURBULENTMAGNETICFIELD_H_ */

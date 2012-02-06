#ifndef TURBULENTMAGNETICFIELD_HPP_
#define TURBULENTMAGNETICFIELD_HPP_

#include "mpc/magneticField/magneticFieldGrid.hpp"

namespace mpc {

/** turbulent field properties
 *  mean and rms strength:
 *  <B> = 0
 *  <B^2> = Brms^2
 *  turbulent range and spectrum:
 *  B(k) ~ k^(-11/3) for 1/8 < k < 1/2
 *  B(k) = 0 otherwise
 *  solenoidity:
 *  B(k) * k = 0
 *  isotropy:
 *  same coherence length in all directions
 */
class TurbulentMagneticFieldGrid: public MagneticFieldGrid {
public:
	TurbulentMagneticFieldGrid(size_t n, double spacing,
			Vector3 origin, double Brms, double lMin, double lMax,
			double powerSpectralIndex, int seed);
protected:
	void initialize();
	double Brms, lMin, lMax, powerSpectralIndex;
	int seed;
};

} // namespace mpc

#endif /* TURBULENTMAGNETICFIELD_HPP_ */

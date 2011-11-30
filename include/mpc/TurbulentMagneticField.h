#ifndef TURBULENTMAGNETICFIELD_H_
#define TURBULENTMAGNETICFIELD_H_

#include "mpc/MagneticField.h"
#include "mpc/Vector3.h"
#include <vector>

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

class TurbulentMagneticField : public MagneticField {
public:
	TurbulentMagneticField(Vector3 origin, size_t n, double spacing,
			double Brms, double spectralIndex, double lMin, double lMax);
	~TurbulentMagneticField();
	Vector3 getField(const Vector3 &position) const;
	void setSeed(unsigned int seed);
	void initialize();

	std::vector<Vector3> field;
	Vector3 origin; // origin of grid
	size_t n; // size of grid
	double spacing; // grid spacing
	double Brms; // root mean square of field strength
	double spectralIndex; // index of turbulent spectrum
	double lMin; // minimum length scale
	double lMax; // maximum length scale
	unsigned int seed;
};

} // namespace mpc

#endif /* TURBULENTMAGNETICFIELD_H_ */

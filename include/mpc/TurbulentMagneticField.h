#ifndef TURBULENTMAGNETICFIELD_H_
#define TURBULENTMAGNETICFIELD_H_

#include "mpc/MagneticField.h"
#include "mpc/Vector3.h"
#include <vector>

namespace mpc {

class TurbulentMagneticField : public MagneticField {
public:
	TurbulentMagneticField(Vector3 origin, size_t N, double spacing,
			double Brms, double spectralIndex, double Lmin, double Lmax);
	~TurbulentMagneticField();
	Vector3 getField(const Vector3 &position) const;
	void initialize();

protected:
	std::vector<Vector3> field;
	Vector3 origin; // origin of grid
	size_t N; // size of grid
	double spacing; // grid spacing
	double Brms; // root mean square of field strength
	double spectralIndex; // index of turbulent spectrum
	double kMin; // minimum length scale
	double kMax; // maximum length scale
	bool periodic; // continue field periodically
};

} // namespace mpc

#endif /* TURBULENTMAGNETICFIELD_H_ */

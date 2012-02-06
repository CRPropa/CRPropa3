#ifndef MAGNETICFIELDGRID_HPP_
#define MAGNETICFIELDGRID_HPP_

#include "mpc/magneticField/magneticField.hpp"
#include "mpc/Vector3.h"
#include <vector>

namespace mpc {

class MagneticFieldGrid: public MagneticField {
public:
	MagneticFieldGrid(size_t n, double spacing, Vector3 origin);
	Vector3 getField(const Vector3 &position) const;

protected:
	std::vector<std::vector<std::vector<Vector3> > > grid;
	size_t n;
	double spacing;
	Vector3 origin;
};

} // namespace mpc

#endif /* MAGNETICFIELDGRID_HPP_ */

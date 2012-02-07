#ifndef MAGNETICFIELDGRID_HPP_
#define MAGNETICFIELDGRID_HPP_

#include "mpc/magneticField/magneticField.hpp"
#include "mpc/Vector3.h"
#include <vector>

namespace mpc {

/**
 @class MagneticFieldGrid
 @brief Cubic, cartesian magnetic field grid with trilinear interpolation.

 This class provides a magnetic field grid. \n
 The grid spacing is constant and equal along all three axes (cartesian). \n
 The grid is of cubic shape. \n
 Magnetic field values are calculated by trilinear interpolation of the surrounding 8 grid points. \n
 The grid is periodically extended.
 */
class MagneticFieldGrid: public MagneticField {
public:
	MagneticFieldGrid(Vector3 origin, size_t n, double spacing);
	Vector3 getField(const Vector3 &position) const;

protected:
	std::vector<std::vector<std::vector<Vector3> > > grid;
	size_t samples;
	double spacing;
	Vector3 origin;
};

} // namespace mpc

#endif /* MAGNETICFIELDGRID_HPP_ */

#ifndef MAGNETICFIELDGRID_H_
#define MAGNETICFIELDGRID_H_

#include "mpc/magneticField/magneticField.h"
#include "mpc/Vector3.h"
#include <vector>

namespace mpc {

/**
 @class MagneticFieldGrid
 @brief Cubic, cartesian magnetic field grid with trilinear interpolation.

 This class provides a magnetic field grid.
 The grid spacing is constant and equal along all three axes (cartesian).
 The grid is of cubic shape.
 Magnetic field values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically extended.
 */
class MagneticFieldGrid: public MagneticField {
public:
	MagneticFieldGrid(Vector3d origin, size_t n, double spacing);
	Vector3d getField(const Vector3d &position) const;
	Vector3d getGridOrigin() const;
	size_t getGridSamples() const;
	double getGridSpacing() const;
	double getGridSize() const;
	virtual void updateSimulationVolume(const Vector3d &origin, double size);

protected:
	std::vector<std::vector<std::vector<Vector3f> > > grid;
	size_t samples;
	double spacing;
	Vector3d origin;
};

} // namespace mpc

#endif /* MAGNETICFIELDGRID_H_ */

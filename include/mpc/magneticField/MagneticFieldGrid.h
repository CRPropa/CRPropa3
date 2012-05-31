#ifndef MPC_MAGNETICFIELDGRID_H_
#define MPC_MAGNETICFIELDGRID_H_

#include "mpc/magneticField/MagneticField.h"
#include <vector>

namespace mpc {

/**
 @class MagneticFieldGrid
 @brief Periodic, cubic, cartesian magnetic field grid with trilinear interpolation.

 This class provides a magnetic field grid.
 The grid spacing is constant and equal along all three axes (cartesian).
 The grid is of cubic shape.
 Magnetic field values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically extended.
 */
class MagneticFieldGrid: public MagneticField {
public:
	MagneticFieldGrid(Vector3d origin, double size, size_t samples);
	Vector3f &get(size_t ix, size_t iy, size_t iz);
	const Vector3f &get(size_t ix, size_t iy, size_t iz) const;
	Vector3d getField(const Vector3d &position) const;
	Vector3d getGridOrigin() const;
	size_t getGridSamples() const;
	double getGridSpacing() const;
	double getGridSize() const;
	virtual void updateSimulationVolume(const Vector3d &origin, double size);

protected:
	/** Grid of magnetic field vectors.
	 *	Since the grid is periodic the field at opposite borders is identical and does not need to be sampled twice.
	 *	The field is thus sampled at 0, size / samples, ... size * (samples - 1) / samples. */
	std::vector<Vector3f> grid;
	Vector3d origin; /** Origin of the field */
	double size; /** Extension of the field */
	size_t samples; /** Number of grid points per edge */
	double spacing; /** Distance between two grid points (= size / samples) */
};

void periodicClamp(double x, int n, int &lo, int &hi);

} // namespace mpc

#endif /* MPC_MAGNETICFIELDGRID_H_ */

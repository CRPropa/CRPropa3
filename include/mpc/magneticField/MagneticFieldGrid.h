#ifndef MPC_MAGNETICFIELDGRID_H_
#define MPC_MAGNETICFIELDGRID_H_

#include "mpc/magneticField/MagneticField.h"
#include <vector>

namespace mpc {

/**
 @class MagneticFieldGrid
 @brief Magnetic field on a periodic, cartesian grid with trilinear interpolation.

 This class provides a three-dimensional magnetic field grid.
 The grid spacing is constant and equal along all three axes.
 Magnetic field values are calculated by trilinear interpolation of the surrounding 8 grid points. \n
 The grid is periodically extended.
 Thus the field at opposite borders is assumed identical and does not need to be sampled twice.
 The grid sample positions are hence at 0, size / samples, ... size * (samples - 1) / samples and the grid spacing is size / samples.
 */
class MagneticFieldGrid: public MagneticField {
	std::vector<Vector3f> grid; /**< Magnetic field vectors */
	Vector3d origin; /**< Lower left front corner of the grid */
	size_t Nx, Ny, Nz; /**< Number of grid points per edge */
	double spacing; /**< Distance between two grid points (= size / samples) */

public:
	MagneticFieldGrid(Vector3d origin, size_t N, double spacing);
	MagneticFieldGrid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz,
			double spacing);

	void setOrigin(Vector3d origin);
	void setGridSize(size_t Nx, size_t Ny, size_t Nz);
	void setSpacing(double spacing);

	Vector3d getOrigin() const;
	size_t getNx() const;
	size_t getNy() const;
	size_t getNz() const;
	double getSpacing() const;

	Vector3f &get(size_t ix, size_t iy, size_t iz);
	const Vector3f &get(size_t ix, size_t iy, size_t iz) const;

	Vector3d getField(const Vector3d &position) const;
};

/** Lower and upper neighbor in a periodically continued unit grid */
void periodicClamp(double x, int n, int &lo, int &hi);

} // namespace mpc

#endif /* MPC_MAGNETICFIELDGRID_H_ */

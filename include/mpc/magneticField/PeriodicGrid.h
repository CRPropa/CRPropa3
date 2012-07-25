#ifndef MPC_PERIODICGRID_H_
#define MPC_PERIODICGRID_H_

#include "mpc/Vector3.h"
#include <vector>

namespace mpc {

// Lower and upper neighbor in a periodically continued unit grid
inline void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ((int(floor(x)) % n) + n) % n;
	hi = (lo + 1) % n;
}

/**
 @class PeriodicGrid
 @brief Template class for fields on a periodic grid with trilinear interpolation

 The grid spacing is constant and equal along all three axes.
 Values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically extended.
 The grid sample positions are at 0, size/N, ... (N-1) * size/N.
 */
template<typename T>
class PeriodicGrid : public Referenced {
	std::vector<T> grid;
	Vector3d origin;
	size_t Nx, Ny, Nz;
	double spacing;

public:
	PeriodicGrid(Vector3d origin, size_t N, double spacing) {
		setOrigin(origin);
		setGridSize(N, N, N);
		setSpacing(spacing);
	}

	PeriodicGrid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz,
			double spacing) {
		setOrigin(origin);
		setGridSize(Nx, Ny, Nz);
		setSpacing(spacing);
	}

	void setOrigin(Vector3d origin) {
		this->origin = origin;
	}

	void setGridSize(size_t Nx, size_t Ny, size_t Nz) {
		this->Nx = Nx;
		this->Ny = Ny;
		this->Nz = Nz;
		grid.resize(Nx * Ny * Nz);
	}

	void setSpacing(double spacing) {
		this->spacing = spacing;
	}

	Vector3d getOrigin() const {
		return origin;
	}
	size_t getNx() const {
		return Nx;
	}

	size_t getNy() const {
		return Ny;
	}

	size_t getNz() const {
		return Nz;
	}

	double getSpacing() const {
		return spacing;
	}

	T &get(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Ny + iz];
	}

	const T &get(size_t ix, size_t iy, size_t iz) const {
		return grid[ix * Ny * Nz + iy * Ny + iz];
	}

	T interpolate(const Vector3d &position) const {
		// position on a unit grid
		Vector3d r = (position - origin) / spacing;

		// indices of lower and upper neighbors
		int ix, iX, iy, iY, iz, iZ;
		periodicClamp(r.x, Nx, ix, iX);
		periodicClamp(r.y, Ny, iy, iY);
		periodicClamp(r.z, Nz, iz, iZ);

		// linear fraction to lower and upper neighbors
		double fx = r.x - floor(r.x);
		double fX = 1 - fx;
		double fy = r.y - floor(r.y);
		double fY = 1 - fy;
		double fz = r.z - floor(r.z);
		double fZ = 1 - fz;

		// trilinear interpolation
		// check: http://paulbourke.net/miscellaneous/interpolation/
		T b(0.);
		//V000 (1 - x) (1 - y) (1 - z) +
		b += get(ix, iy, iz) * fX * fY * fZ;
		//V100 x (1 - y) (1 - z) +
		b += get(iX, iy, iz) * fx * fY * fZ;
		//V010 (1 - x) y (1 - z) +
		b += get(ix, iY, iz) * fX * fy * fZ;
		//V001 (1 - x) (1 - y) z +
		b += get(ix, iy, iZ) * fX * fY * fz;
		//V101 x (1 - y) z +
		b += get(iX, iy, iZ) * fx * fY * fz;
		//V011 (1 - x) y z +
		b += get(ix, iY, iZ) * fX * fy * fz;
		//V110 x y (1 - z) +
		b += get(iX, iY, iz) * fx * fy * fZ;
		//V111 x y z
		b += get(iX, iY, iZ) * fx * fy * fz;

		return b;
	}
};

typedef PeriodicGrid<Vector3f> VectorGrid;
typedef PeriodicGrid<float> ScalarGrid;

} // namespace mpc

#endif /* MPC_PERIODICGRID_H_ */

#ifndef CRPROPA_GRID_H
#define CRPROPA_GRID_H

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include <vector>

namespace crpropa {

/** Lower and upper neighbor in a periodically continued unit grid */
inline void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ((int(floor(x)) % n) + n) % n;
	hi = (lo + 1) % n;
}

/** Lower and upper neighbor in a reflectively repeated unit grid */
inline void reflectiveClamp(double x, int n, int &lo, int &hi) {
	while ((x < 0) or (x > n))
		x = 2 * n * (x > n) - x;
	lo = floor(x);
	hi = lo + (lo < n-1);
}

/** Symmetrical round */
inline double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

/**
 * \addtogroup Core
 * @{
 */
/**
 @class Grid
 @brief Template class for fields on a periodic grid with trilinear interpolation

 The grid spacing is constant with diffrent resulution along all three axes.
 Values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically (default) or reflectively extended.
 The grid sample positions are at 1/2 * size/N, 3/2 * size/N ... (2N-1)/2 * size/N.
 */
template<typename T>
class Grid: public Referenced {
	std::vector<T> grid;
	size_t Nx, Ny, Nz; /**< Number of grid points */
	Vector3d origin; /**< Origin of the volume that is represented by the grid. */
	Vector3d gridOrigin; /**< Grid origin */
	Vector3d spacing; /**< Distance between grid points, determines the extension of the grid */
	bool reflective; /**< If set to true, the grid is repeated reflectively instead of periodically */

public:
	/** Constructor for cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	N		Number of grid points in one direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t N, double spacing) {
		setOrigin(origin);
		setGridSize(N, N, N);
		setSpacing(Vector3d(spacing));
		setReflective(false);
	}

	/** Constructor for non-cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) {
		setOrigin(origin);
		setGridSize(Nx, Ny, Nz);
		setSpacing(Vector3d(spacing));
		setReflective(false);
	}
	
	/** Constructor for non-cubic grid with spacing vector
	@param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing vector between grid points
	*/
	 Grid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, Vector3d spacing) {
	 	setOrigin(origin);
	 	setGridSize(Nx, Ny, Nz);
	 	setSpacing(spacing);
	 	setReflective(false);
	 } 

	void setOrigin(Vector3d origin) {
		this->origin = origin;
		this->gridOrigin = origin + spacing/2;
	}

	/** Resize grid, also enlarges the volume as the spacing stays constant */
	void setGridSize(size_t Nx, size_t Ny, size_t Nz) {
		this->Nx = Nx;
		this->Ny = Ny;
		this->Nz = Nz;
		grid.resize(Nx * Ny * Nz);
		setOrigin(origin);
	}

	void setSpacing(Vector3d spacing) {
		this->spacing = spacing;
		setOrigin(origin);
	}

	void setReflective(bool b) {
		reflective = b;
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

	Vector3d getSpacing() const {
		return spacing;
	}

	bool isReflective() const {
		return reflective;
	}

	/** Inspector & Mutator */
	T &get(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	/** Inspector */
	const T &get(size_t ix, size_t iy, size_t iz) const {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	T getValue(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	void setValue(size_t ix, size_t iy, size_t iz, T value) {
		grid[ix * Ny * Nz + iy * Nz + iz] = value;
	}

	/** Return a reference to the grid values */
	std::vector<T> &getGrid() {
		return grid;
	}

	/** Position of the grid point of a given index */
	Vector3d positionFromIndex(int index) const {
		int ix = index / (Ny * Nz);
		int iy = (index / Nz) % Ny;
		int iz = index % Nz;
		return Vector3d(ix, iy, iz) * spacing + gridOrigin;
	}

	/** Value of a grid point that is closest to a given position */
	T closestValue(const Vector3d &position) const {
		Vector3d r = (position - gridOrigin) / spacing;
		int ix = round(r.x);
		int iy = round(r.y);
		int iz = round(r.z);
		if (reflective) {
			while ((ix < 0) or (ix > Nx))
				ix = 2 * Nx * (ix > Nx) - ix;
			while ((iy < 0) or (iy > Ny))
				iy = 2 * Ny * (iy > Ny) - iy;
			while ((iz < 0) or (iz > Nz))
				iz = 2 * Nz * (iz > Nz) - iz;
		} else {
			ix = ((ix % Nx) + Nx) % Nx;
			iy = ((iy % Ny) + Ny) % Ny;
			iz = ((iz % Nz) + Nz) % Nz;
		}
		return get(ix, iy, iz);
	}

	/** Interpolate the grid at a given position */
	T interpolate(const Vector3d &position) const {
		// position on a unit grid
		Vector3d r = (position - gridOrigin) / spacing;

		// indices of lower and upper neighbors
		int ix, iX, iy, iY, iz, iZ;
		if (reflective) {
			reflectiveClamp(r.x, Nx, ix, iX);
			reflectiveClamp(r.y, Ny, iy, iY);
			reflectiveClamp(r.z, Nz, iz, iZ);
		} else {
			periodicClamp(r.x, Nx, ix, iX);
			periodicClamp(r.y, Ny, iy, iY);
			periodicClamp(r.z, Nz, iz, iZ);
		}

		// linear fraction to lower and upper neighbors
		double fx = r.x - floor(r.x);
		double fX = 1 - fx;
		double fy = r.y - floor(r.y);
		double fY = 1 - fy;
		double fz = r.z - floor(r.z);
		double fZ = 1 - fz;

		// trilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation)
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

typedef Grid<Vector3f> VectorGrid;
typedef Grid<float> ScalarGrid;
/** @}*/

} // namespace crpropa

#endif // CRPROPA_GRID_H

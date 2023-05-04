#ifndef CRPROPA_GRID_H
#define CRPROPA_GRID_H

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

#include "kiss/string.h"
#include "kiss/logger.h"

#include <vector>
#include <type_traits>
#if HAVE_SIMD
#include <immintrin.h>
#include <smmintrin.h>
#endif // HAVE_SIMD

namespace crpropa {

/** If set to TRILINEAR, use trilinear interpolation (standard)
If set to TRICUBIC, use tricubic interpolation instead of trilinear interpolation
If set to NEAREST_NEIGHBOUR , use nearest neighbour interpolation instead of trilinear interpolation */
enum interpolationType {
  TRILINEAR = 0,
  TRICUBIC,
  NEAREST_NEIGHBOUR
};

/** Lower and upper neighbour in a periodically continued unit grid */
inline void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ((int(floor(x)) % (n)) + (n)) % (n);
	hi = (lo + 1) % (n);
}

/** grid index in a reflective continued unit grid */
inline int reflectiveBoundary(int index, int n) {
	while ((index < -0.5) or (index > (n-0.5)))
		index = 2 * n * (index > (n-0.5)) - index-1;
	return index;
}

/** grid index in a periodically continued unit grid */
inline int periodicBoundary(int index, int n) {
	return ((index % (n)) + (n)) % (n);
}

/** Lower and upper neighbour in a reflectively repeated unit grid */
inline void reflectiveClamp(double x, int n, int &lo, int &hi, double &res) {
	while ((x < -0.5) or (x > (n-0.5)))
		x = 2 * n * (x > (n-0.5)) -x-1;
	res = x;
	lo = floor(x);
	hi = lo + (lo < n-1);
	if (x<0) {
		lo=0; 
		hi=0;
	}
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
 @class GridProperties
 @brief Combines parameters that uniquely define Grid class
 */
class GridProperties: public Referenced {
public:
	size_t Nx, Ny, Nz; 	// Number of grid points
	Vector3d origin;  	// Position of the lower left front corner of the volume
	Vector3d spacing; 	// Spacing vector between gridpoints
	bool reflective;	// using reflective repetition of the grid instead of periodic
	interpolationType ipol;	// Interpolation type used between grid points
	bool clipVolume;	// Set grid values to 0 outside the volume if true

	/** Constructor for cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	N		Number of grid points in one direction
	 @param spacing	Spacing between grid points
	 */
	GridProperties(Vector3d origin, size_t N, double spacing) :
		origin(origin), Nx(N), Ny(N), Nz(N), spacing(Vector3d(spacing)), reflective(false), ipol(TRILINEAR), clipVolume(false) {
	}

	/** Constructor for non-cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing between grid points
	 */
	GridProperties(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) :
		origin(origin), Nx(Nx), Ny(Ny), Nz(Nz), spacing(Vector3d(spacing)), reflective(false), ipol(TRILINEAR), clipVolume(false) {
	}

	/** Constructor for non-cubic grid with spacing vector
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing vector between grid points
	*/
	GridProperties(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, Vector3d spacing) :
		origin(origin), Nx(Nx), Ny(Ny), Nz(Nz), spacing(spacing), reflective(false), ipol(TRILINEAR), clipVolume(false) {
	}
	
	virtual ~GridProperties() {
	}
	
	/** If True, the repetition of the grid is refletive instead of periodic. */
	void setReflective(bool b) {
		reflective = b;
	}

	/** set the type of interpolation between grid points.
	 * @param i: interpolationType (TRILINEAR, TRICUBIC, NEAREST_NEIGHBOUR) */
	void setInterpolationType(interpolationType i) {
		ipol = i;
	}

	/** If True, the grid is set to zero outside of the volume. */
	void setClipVolume(bool b) {
		clipVolume = b;
	}
};

/**
 @class Grid
 @brief Template class for fields on a periodic grid with trilinear interpolation

 The grid spacing is constant with diffrent resolution along all three axes.
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
	bool clipVolume; /**< If set to true, all values outside of the grid will be 0*/
	bool reflective; /**< If set to true, the grid is repeated reflectively instead of periodically */
	interpolationType ipolType; /**< Type of interpolation between the grid points */

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
		setClipVolume(false);
		setInterpolationType(TRILINEAR);
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
		setClipVolume(false);
		setInterpolationType(TRILINEAR);
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
		setClipVolume(false);
		setInterpolationType(TRILINEAR);
	}

	/** Constructor for GridProperties
	 @param p	GridProperties instance
     */
	Grid(const GridProperties &p) :
		origin(p.origin), spacing(p.spacing), reflective(p.reflective), ipolType(p.ipol) {
		setGridSize(p.Nx, p.Ny, p.Nz);
		setClipVolume(p.clipVolume);
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

	// If set to true, all values outside of the grid will be 0.
	void setClipVolume(bool b) {
		clipVolume = b;
	}

	/** Change the interpolation type to the routine specified by the user. Check if this routine is
		contained in the enum interpolationType and thus supported by CRPropa.*/
	void setInterpolationType(interpolationType ipolType) {
		if (ipolType == TRILINEAR || ipolType == TRICUBIC || ipolType == NEAREST_NEIGHBOUR) {
			this->ipolType = ipolType;
			if ((ipolType == TRICUBIC) && (std::is_same<T, Vector3d>::value)) {
				KISS_LOG_WARNING << "Tricubic interpolation on Grid3d works only with float-precision, doubles will be downcasted";
		}
		} else {
			throw std::runtime_error("InterpolationType: unknown interpolation type");
		}
	}

	/** returns the position of the lower left front corner of the volume */
	Vector3d getOrigin() const {
		return origin;
	}

	bool getClipVolume() const {
		return clipVolume;
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

	/** Calculates the total size of the grid in bytes */
	size_t getSizeOf() const {
		return sizeof(grid) + (sizeof(grid[0]) * grid.size());
	}

	Vector3d getSpacing() const {
		return spacing;
	}

	bool isReflective() const {
		return reflective;
	}

	/** Choose the interpolation algorithm based on the set interpolation type.
	  By default this it the trilinear interpolation. The user can change the
	  routine with the setInterpolationType function.*/
	T interpolate(const Vector3d &position) {
		// check for volume
		if (clipVolume) {
			Vector3d edge = origin + Vector3d(Nx, Ny, Nz) * spacing;
			bool isInVolume = (position.x >= origin.x) && (position.x <= edge.x);
			isInVolume &= (position.y >= origin.y) && (position.y <= edge.y);
			isInVolume &= (position.z >= origin.z) && (position.z <= edge.z);
			if (!isInVolume) 
				return T(0.);
		} 

		if (ipolType == TRICUBIC)
			return tricubicInterpolate(T(), position);
		else if (ipolType == NEAREST_NEIGHBOUR)
			return closestValue(position);
		else
			return trilinearInterpolate(position);
	}

	/** Inspector & Mutator */
	T &get(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	/** Inspector */
	const T &get(size_t ix, size_t iy, size_t iz) const {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	const T &periodicGet(size_t ix, size_t iy, size_t iz) const {
		ix = periodicBoundary(ix, Nx);
		iy = periodicBoundary(iy, Ny);
		iz = periodicBoundary(iz, Nz);
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	const T &reflectiveGet(size_t ix, size_t iy, size_t iz) const {
		ix = reflectiveBoundary(ix, Nx);
		iy = reflectiveBoundary(iy, Ny);
		iz = reflectiveBoundary(iz, Nz);
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

	/** Value of a grid point that is closest to a given position / nearest neighbour interpolation */
	T closestValue(const Vector3d &position) const {
		Vector3d r = (position - gridOrigin) / spacing;
		int ix, iy, iz;
		if (reflective) {
			ix = round(r.x);
			iy = round(r.y);
			iz = round(r.z);
			while ((ix < -0.5) or (ix > (Nx-0.5)))
				ix = 2 * Nx * (ix > (Nx-0.5)) - ix-1;
			while ((iy < -0.5) or (iy > (Ny-0.5)))
				iy = 2 * Ny * (iy > (Ny-0.5)) - iy-1;
			while ((iz < -0.5) or (iz > (Nz-0.5)))
				iz = 2 * Nz * (iz > (Nz-0.5)) - iz-1;
		} else {
			ix = round(fmod(r.x, Nx));
			iy = round(fmod(r.y, Ny));
			iz = round(fmod(r.z, Nz));
			ix = (ix + Nx * (ix < 0)) % Nx;
			iy = (iy + Ny * (iy < 0)) % Ny;
			iz = (iz + Nz * (iz < 0)) % Nz;
		}
		return get(ix, iy, iz);
	}

private:
	#ifdef HAVE_SIMD
	__m128 simdperiodicGet(size_t ix, size_t iy, size_t iz) const {
		ix = periodicBoundary(ix, Nx);
		iy = periodicBoundary(iy, Ny);
		iz = periodicBoundary(iz, Nz);
		return convertVector3fToSimd(grid[ix * Ny * Nz + iy * Nz + iz]);
	}

	__m128 simdreflectiveGet(size_t ix, size_t iy, size_t iz) const {
		ix = reflectiveBoundary(ix, Nx);
		iy = reflectiveBoundary(iy, Ny);
		iz = reflectiveBoundary(iz, Nz);
		return convertVector3fToSimd(grid[ix * Ny * Nz + iy * Nz + iz]);
	}

	__m128 convertVector3fToSimd(const Vector3f v) const {
		__m128 simdVar = _mm_set_ps(0,v.z,v.y,v.x);
		return simdVar;
	}
	
	Vector3f convertSimdToVector3f(__m128 res) const {
		float vec[4];	
		_mm_store_ps(&vec[0], res);
		Vector3f result = Vector3f(vec[0], vec[1], vec[2]);
		return result;
	}

	/** Vectorized cubic Interpolator in 1D */
	__m128 CubicInterpolate(__m128 p0,__m128 p1,__m128 p2,__m128 p3,double position) const {
		__m128 c1 = _mm_set1_ps (1/2.);
		__m128 c2 = _mm_set1_ps (3/2.);
		__m128 c3 = _mm_set1_ps (2.);
		__m128 c4 = _mm_set1_ps (5/2.);

		__m128 pos  = _mm_set1_ps (position);
		__m128 pos2 = _mm_set1_ps (position*position);
		__m128 pos3 = _mm_set1_ps (position*position*position);

		/** SIMD optimized routine to calculate 'res = ((-0.5*p0+3/2.*p1-3/2.*p2+0.5*p3)*pos*pos*pos+(p0-5/2.*p1+p2*2-0.5*p3)*pos*pos+(-0.5*p0+0.5*p2)*pos+p1);'
			 where terms are used as:
			term = (-0.5*p0+0.5*p2)*pos
			term2 = (p0-5/2.*p1+p2*2-0.5*p3)*pos*pos;
			term3 = (-0.5*p0+3/2.*p1-3/2.*p2+0.5*p3)*pos*pos*pos;  */
		__m128 term = _mm_mul_ps(_mm_sub_ps(_mm_mul_ps(c1,p2),_mm_mul_ps(c1,p0)),pos);
		__m128 term2 = _mm_mul_ps(_mm_sub_ps(_mm_add_ps(p0,_mm_mul_ps(c3,p2)),_mm_add_ps(_mm_mul_ps(c4,p1),_mm_mul_ps(c1,p3))),pos2);
		__m128 term3 = _mm_mul_ps(_mm_sub_ps(_mm_add_ps(_mm_mul_ps(c2,p1),_mm_mul_ps(c1,p3)),_mm_add_ps(_mm_mul_ps(c1,p0),_mm_mul_ps(c2,p2))),pos3);
		__m128 res = _mm_add_ps(_mm_add_ps(_mm_add_ps(term3,term2),term),p1);
		return res;
	}
	#endif // HAVE_SIMD
	/** Interpolate the grid tricubic at a given position (see https://www.paulinternet.nl/?page=bicubic, http://graphics.cs.cmu.edu/nsp/course/15-462/Fall04/assts/catmullRom.pdf) */
	Vector3f tricubicInterpolate(Vector3f, const Vector3d &position) const {
		#ifdef HAVE_SIMD
		// position on a unit grid
		Vector3d r = (position - gridOrigin) / spacing;

		int iX0, iY0, iZ0;
		iX0 = floor(r.x);
		iY0 = floor(r.y);
		iZ0 = floor(r.z);

		double fX, fY, fZ;
		fX = r.x - iX0;
		fY = r.y - iY0;
		fZ = r.z - iZ0;

		int nrCubicInterpolations = 4;
		__m128 interpolateVaryX[nrCubicInterpolations];
		__m128 interpolateVaryY[nrCubicInterpolations];
		__m128 interpolateVaryZ[nrCubicInterpolations];
		/** Perform 1D interpolations while iterating in each for loop over the index of another direction */
		for (int iLoopX = -1; iLoopX < nrCubicInterpolations-1; iLoopX++) {
			for (int iLoopY = -1; iLoopY < nrCubicInterpolations-1; iLoopY++) {
				for (int iLoopZ = -1; iLoopZ < nrCubicInterpolations-1; iLoopZ++) {
					if (reflective)
						interpolateVaryZ[iLoopZ+1] = simdreflectiveGet(iX0+iLoopX, iY0+iLoopY, iZ0+iLoopZ);
					else 
						interpolateVaryZ[iLoopZ+1] = simdperiodicGet(iX0+iLoopX, iY0+iLoopY, iZ0+iLoopZ);
				}
				interpolateVaryY[iLoopY+1] = CubicInterpolate(interpolateVaryZ[0], interpolateVaryZ[1], interpolateVaryZ[2], interpolateVaryZ[3], fZ);
			}
			interpolateVaryX[iLoopX+1] = CubicInterpolate(interpolateVaryY[0], interpolateVaryY[1], interpolateVaryY[2], interpolateVaryY[3], fY);
		}
		__m128 result = CubicInterpolate(interpolateVaryX[0], interpolateVaryX[1], interpolateVaryX[2], interpolateVaryX[3], fX);
		return convertSimdToVector3f(result);
		#else // HAVE_SIMD
		throw std::runtime_error( "Tried to use tricubic Interpolation without SIMD_EXTENSION. SIMD Optimization is necessary for tricubic interpolation of vector grids.\n");
		#endif // HAVE_SIMD	
	}

	/** Vectorized cubic Interpolator in 1D that returns a scalar (see https://www.paulinternet.nl/?page=bicubic, http://graphics.cs.cmu.edu/nsp/course/15-462/Fall04/assts/catmullRom.pdf) */
	double CubicInterpolateScalar(double p0,double p1,double p2,double p3,double pos) const {
		return((-0.5*p0+3/2.*p1-3/2.*p2+0.5*p3)*pos*pos*pos+(p0-5/2.*p1+p2*2-0.5*p3)*pos*pos+(-0.5*p0+0.5*p2)*pos+p1);
	}

  /** Interpolate the grid tricubic at a given position (see https://www.paulinternet.nl/?page=bicubic, http://graphics.cs.cmu.edu/nsp/course/15-462/Fall04/assts/catmullRom.pdf) */
	double tricubicInterpolate(double, const Vector3d &position) const {
		/** position on a unit grid */
		Vector3d r = (position - gridOrigin) / spacing;

		int iX0, iY0, iZ0;
		iX0 = floor(r.x);
		iY0 = floor(r.y);
		iZ0 = floor(r.z);

		double fX, fY, fZ;
		fX = r.x - iX0;
		fY = r.y - iY0;
		fZ = r.z - iZ0;

		int nrCubicInterpolations = 4;
		double interpolateVaryX[nrCubicInterpolations];
		double interpolateVaryY[nrCubicInterpolations];
		double interpolateVaryZ[nrCubicInterpolations];
		/** Perform 1D interpolations while iterating in each for loop over the index of another direction */
		for (int iLoopX = -1; iLoopX < nrCubicInterpolations-1; iLoopX++) {
			for (int iLoopY = -1; iLoopY < nrCubicInterpolations-1; iLoopY++) {
				for (int iLoopZ = -1; iLoopZ < nrCubicInterpolations-1; iLoopZ++) {
					if (reflective)
						interpolateVaryZ[iLoopZ+1] = reflectiveGet(iX0+iLoopX, iY0+iLoopY, iZ0+iLoopZ);
					else
						interpolateVaryZ[iLoopZ+1] = periodicGet(iX0+iLoopX, iY0+iLoopY, iZ0+iLoopZ);
				}
				interpolateVaryY[iLoopY+1] = CubicInterpolateScalar(interpolateVaryZ[0], interpolateVaryZ[1], interpolateVaryZ[2], interpolateVaryZ[3], fZ);
			}
			interpolateVaryX[iLoopX+1] = CubicInterpolateScalar(interpolateVaryY[0], interpolateVaryY[1], interpolateVaryY[2], interpolateVaryY[3], fY);
		}
		double result = CubicInterpolateScalar(interpolateVaryX[0], interpolateVaryX[1], interpolateVaryX[2], interpolateVaryX[3], fX);
		return result;
	}

	/** Interpolate the grid trilinear at a given position */
	T trilinearInterpolate(const Vector3d &position) const {
		/** position on a unit grid */
		Vector3d r = (position - gridOrigin) / spacing;

		/** indices of lower (0) and upper (1) neighbours. The neighbours span a grid
		  with the origin at [iX0, iY0, iZ0] and the most distant corner [iX1, iY1, iZ1]. */
		int iX0, iX1, iY0, iY1, iZ0, iZ1;
		double resX, resY, resZ, fX0, fY0, fZ0;

		if (reflective) {
			reflectiveClamp(r.x, Nx, iX0, iX1, resX);
			reflectiveClamp(r.y, Ny, iY0, iY1, resY);
			reflectiveClamp(r.z, Nz, iZ0, iZ1, resZ);
			fX0 = resX - floor(resX);
			fY0 = resY - floor(resY);
			fZ0 = resZ - floor(resZ);
		} else {
			periodicClamp(r.x, Nx, iX0, iX1);
			periodicClamp(r.y, Ny, iY0, iY1);
			periodicClamp(r.z, Nz, iZ0, iZ1);
			fX0 = r.x - floor(r.x);
			fY0 = r.y - floor(r.y);
			fZ0 = r.z - floor(r.z);
		}

		/** linear fraction to upper neighbours based on lower neighbours calculated above */
		double fX1 = 1 - fX0;
		double fY1 = 1 - fY0;
		double fZ1 = 1 - fZ0;

		/** trilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation) */
		T b(0.);
		b += get(iX0, iY0, iZ0) * fX1 * fY1 * fZ1;
		b += get(iX1, iY0, iZ0) * fX0 * fY1 * fZ1;
		b += get(iX0, iY1, iZ0) * fX1 * fY0 * fZ1;
		b += get(iX0, iY0, iZ1) * fX1 * fY1 * fZ0;
		b += get(iX1, iY0, iZ1) * fX0 * fY1 * fZ0;
		b += get(iX0, iY1, iZ1) * fX1 * fY0 * fZ0;
		b += get(iX1, iY1, iZ0) * fX0 * fY0 * fZ1;
		b += get(iX1, iY1, iZ1) * fX0 * fY0 * fZ0;

		return b;
	}

}; // class Grid

typedef Grid<double> Grid1d;
typedef Grid<float> Grid1f;
typedef Grid<Vector3f> Grid3f;
typedef Grid<Vector3d> Grid3d;

/** @}*/

} // namespace crpropa

#endif // CRPROPA_GRID_H

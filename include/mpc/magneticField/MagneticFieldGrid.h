#ifndef MPC_MAGNETICFIELDGRID_H_
#define MPC_MAGNETICFIELDGRID_H_

#include "mpc/magneticField/MagneticField.h"

#include <vector>

namespace mpc {

/**
 @class MagneticFieldGrid
 @brief Periodic, cubic, cartesian magnetic field grid with trilinear interpolation.

 This class provides a magnetic field grid.
 The grid is of cubic shape and the spacing is constant and equal along all three axes.
 Magnetic field values are calculated by trilinear interpolation of the surrounding 8 grid points. \n
 The grid is periodically extended.
 Thus the field at opposite borders is assumed identical and does not need to be sampled twice.
 The grid sample positions are hence at 0, size / samples, ... size * (samples - 1) / samples and the grid spacing is size / samples.
 */
class MagneticFieldGrid: public MagneticField {
public:
	/**
	 * Constructor
	 * @param origin	Lower left front corner of the grid
	 * @param size 	Physical extension of the field grid
	 * @param samples	Number of grid samples per edge
	 */
	MagneticFieldGrid(Vector3d origin, double size, size_t samples);

	/**
	 * Multiply the field by a factor.
	 * @param norm	Normalization factor
	 */
	void normalize(double norm);

	/**
	 * Load the field from a binary file.
	 * The field is stored single precision numbers with the field components in xyz order and the grid z-index changing the fastest.
	 */
	void load(std::string filename);

	/**
	 * Dump the field to a binary file.
	 * The field is stored single precision numbers with the field components in xyz order and the grid z-index changing the fastest.
	 */
	void dump(std::string filename) const;

	/**
	 * Modulate the magnetic field with a density field from a binary file.
	 * The field should be (re)normalized afterwards.
	 * The file has to contain the density values in float precision on a grid with N + 1 samples per edge, the last of which is disregarded due to periodicity.
	 * @param filename	Path to the density grid file
	 * @param exp		Exponent for the modulation
	 */
	void modulateWithDensityField(std::string filename, double exp = 2. / 3.);

	/** Return a reference to the grid point (ix, iy, iz). */
	Vector3f &get(size_t ix, size_t iy, size_t iz);
	const Vector3f &get(size_t ix, size_t iy, size_t iz) const;

	/** Return the field vector at arbitrary position by trilinear interpolation. */
	Vector3d getField(const Vector3d &position) const;

	/** Return the RMS field strength */
	double getRMSFieldStrength() const;

	/**
	 * Return the RMS field strength inside a spherical region
	 * @param center	Center of the sphere
	 * @param radius	Radius of the sphere
	 */
	double getRMSFieldStrengthInSphere(Vector3d center, double radius) const;

	Vector3d getGridOrigin() const;
	size_t getGridSamples() const;
	double getGridSpacing() const;
	double getGridSize() const;

	/** Update the fields origin and size. */
	virtual void updateSimulationVolume(const Vector3d &origin, double size);


protected:
	std::vector<Vector3f> grid; /**< Grid of magnetic field vectors */
	Vector3d origin; /**< Origin of the field */
	double size; /**< Extension of the field */
	size_t samples; /**< Number of grid points per edge */

	double spacing; /**< Distance between two grid points (= size / samples) */
};

/** Lower and upper neighbor in a periodically continued unit grid */
void periodicClamp(double x, int n, int &lo, int &hi);

} // namespace mpc

#endif /* MPC_MAGNETICFIELDGRID_H_ */

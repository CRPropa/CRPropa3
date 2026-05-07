#ifndef CRPROPA_GRIDTOOLS_H
#define CRPROPA_GRIDTOOLS_H

#include "crpropa/Grid.h"
#include "crpropa/magneticField/MagneticField.h"
#include <string>
#include <array>

#ifdef CRPROPA_HAVE_FFTW3F
#include "fftw3.h"
#endif

/**
 @file
 @brief Grid related functions: load, dump, save, retrieve grid properties ...

 This file contains a number of functions related to scalar and vector grids (Grid.h).

 Dump/load functions are available for saving/loading grids to/from and binary and plain text files.
 In the files the grid points are stored from (0, 0, 0) to (Nx, Ny, Nz) with the z-index changing the fastest.
 Vector components are stored per grid point in xyz-order.
 In case of plain-text files the vector components are separated by a blank or tab and grid points are stored one per line.
 All functions offer a conversion factor that is multiplied to all values.
 */

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/** Evaluate the mean vector of all grid points.
 @param grid		a vector grid (Grid3f)
 @returns The mean of all grid points along each axis
 */
Vector3f meanFieldVector(ref_ptr<Grid3f> grid);

/** Evaluate the mean of all grid points. 
 @param grid		a scalar grid (Grid1f)
 @returns The mean of all grid points
 */
double meanFieldStrength(ref_ptr<Grid1f> grid);
/** Evaluate the mean of all grid points.
 @param grid		a vector grid (Grid3f)
 @returns The mean of all grid points
 */
double meanFieldStrength(ref_ptr<Grid3f> grid);

/** Evaluate the RMS of all grid points.
 @param grid		a scalar grid (Grid1f)
 @returns The total RMS of all grid points.
 */
double rmsFieldStrength(ref_ptr<Grid1f> grid);
/** Evaluate the RMS of all grid points.
 @param grid		a vector grid (Grid3f)
 @returns The total RMS of all grid points.
 */
double rmsFieldStrength(ref_ptr<Grid3f> grid);
/** Evaluate the RMS of all grid points per axis. 
 @param grid		a vector grid (Grid3f)
 @returns An array of length 3 with the RMS field along each axis.
 */
std::array<float, 3> rmsFieldStrengthPerAxis(ref_ptr<Grid3f> grid);

/** Multiply all grid values by a given factor.
 @param grid		a scalar grid (Grid1f)
 @param a			scaling factor that will be used to multiply all points in grid
 */
void scaleGrid(ref_ptr<Grid1f> grid, double a);
/** Multiply all grid values by a given factor.
 @param grid		a vector grid (Grid3f)
 @param a			scaling factor that will be used to multiply all points in grid
 */
void scaleGrid(ref_ptr<Grid3f> grid, double a);

/** Fill vector grid from provided magnetic field.
 @param grid		a vector grid (Grid3f)
 @param field		the magnetic field 
 */
void fromMagneticField(ref_ptr<Grid3f> grid, ref_ptr<MagneticField> field);

/** Fill scalar grid from provided magnetic field.
 @param grid		a scalar grid (Grid1f)
 @param field		the magnetic field
 */
void fromMagneticFieldStrength(ref_ptr<Grid1f> grid, ref_ptr<MagneticField> field);

/** Load a Grid3f from a binary file with single precision.
 @param grid		a vector grid (Grid3f)
 @param filename	name of input file
 @param conversion	multiply every point in grid by a conversion factor
 */
void loadGrid(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Load a Grid1f from a binary file with single precision.
 @param grid		a scalar grid (Grid1f)
 @param filename	name of input file
 @param conversion	multiply every point in grid by a conversion factor
 */
void loadGrid(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid3f to a binary file.
 @param grid		a vector grid (Grid3f)
 @param filename	name of input file
 @param conversion	multiply every point in grid by a conversion factor
 */
void dumpGrid(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid1f to a binary file with single precision.
 @param grid		a scalar grid (Grid1f)
 @param filename	name of input file
 @param conversion	multiply every point in grid by a conversion factor
 */
void dumpGrid(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** Load a Grid3f from a plain text file based on the gridproperties stored in the header
 @param filename	name of the input file
 @param conversion	multiply every point in a grid by a conversion factor
*/
ref_ptr<Grid3f> loadGrid3fFromTxt(std::string filename, double conversion = 1);

/** Load a Grid3f grid from a plain text file.
 @param grid		a vector grid (Grid3f) to which the points will be loaded
 @param filename	name of input file
 @param conversion	multiply every point in grid by a conversion factor
 */
void loadGridFromTxt(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Load a Grid1f from a plain text file based on the gridproperties stored in the header
 @param filename	name of the input file
 @param conversion	multiply every point in a grid by a conversion factor
*/
ref_ptr<Grid1f> loadGrid1fFromTxt(std::string filename, double conversion = 1);

/** Load a Grid1f from a plain text file
 @param grid		a scalar grid (Grid1f) to which the points will be loaded
 @param filename	name of input file
 @param conversion	multiply every point in grid by a conversion factor
 */
void loadGridFromTxt(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid3f to a plain text file.
 @param grid		a vector grid (Grid3f)
 @param filename	name of output file
 @param conversion	multiply every point in grid by a conversion factor
 @param storeProperties	if true the grid properties are stored as a comment
 */
void dumpGridToTxt(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1, bool storeProperties = false);

/** Dump a Grid1f to a plain text file. 
 @param grid		a scalar grid (Grid1f)
 @param filename	name of output file
 @param conversion	multiply every point in grid by a conversion factor
 @param storeProperties	if true the grid properties are stored as a comment
 */
void dumpGridToTxt(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1, bool storeProperties = false);

#ifdef CRPROPA_HAVE_FFTW3F
/**
 Calculate the omnidirectional power spectrum E(k) for a given turbulent field
 @param grid	a three-dimensional grid
 @returns Returns a vector of pairs (k_i, E(k_i)).
*/
std::vector<std::pair<int, float>> gridPowerSpectrum(ref_ptr<Grid3f> grid);
#endif // CRPROPA_HAVE_FFTW3F

/** @}*/
} // namespace crpropa

#endif // CRPROPA_GRIDTOOLS_H

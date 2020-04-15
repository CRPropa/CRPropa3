#ifndef CRPROPA_GRIDTOOLS_H
#define CRPROPA_GRIDTOOLS_H

#include "crpropa/Grid.h"
#include "crpropa/magneticField/MagneticField.h"
#include <string>

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

/** Evaluate the mean vector of all grid points */
Vector3f meanFieldVector(ref_ptr<Grid3f> grid);

/** Evaluate the mean of all grid points */
double meanFieldStrength(ref_ptr<Grid1f> grid);
/** Evaluate the mean of all grid points */
double meanFieldStrength(ref_ptr<Grid3f> grid);

/** Evaluate the RMS of all grid points */
double rmsFieldStrength(ref_ptr<Grid1f> grid);
/** Evaluate the RMS of all grid points */
double rmsFieldStrength(ref_ptr<Grid3f> grid);

/** Multiply all grid values by a given factor */
void scaleGrid(ref_ptr<Grid1f> grid, double a);
/** Multiply all grid values by a given factor */
void scaleGrid(ref_ptr<Grid3f> grid, double a);

/** Fill vector grid from provided magnetic field */
void fromMagneticField(ref_ptr<Grid3f> grid, ref_ptr<MagneticField> field);

/** Fill scalar grid from provided magnetic field */
void fromMagneticFieldStrength(ref_ptr<Grid1f> grid, ref_ptr<MagneticField> field);

/** Load a Grid3f from a binary file with single precision */
void loadGrid(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Load a Grid1f from a binary file with single precision */
void loadGrid(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid3f to a binary file */
void dumpGrid(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid1f to a binary file with single precision */
void dumpGrid(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** Load a Grid3f grid from a plain text file */
void loadGridFromTxt(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Load a Grid1f from a plain text file */
void loadGridFromTxt(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid3f to a plain text file */
void dumpGridToTxt(ref_ptr<Grid3f> grid, std::string filename,
		double conversion = 1);

/** Dump a Grid1f to a plain text file */
void dumpGridToTxt(ref_ptr<Grid1f> grid, std::string filename,
		double conversion = 1);

/** @}*/
} // namespace crpropa

#endif // CRPROPA_GRIDTOOLS_H

#ifndef CRPROPA_GRIDTOOLS_H
#define CRPROPA_GRIDTOOLS_H

#include "crpropa/Grid.h"
#include "crpropa/magneticField/MagneticField.h"
#include <string>

/**
 @file
 @brief Grid related functions: load, dump, save, create turbulent field ...

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
Vector3f meanFieldVector(ref_ptr<VectorGrid> grid);

/** Evaluate the mean of all grid points */
double meanFieldStrength(ref_ptr<ScalarGrid> grid);
/** Evaluate the mean of all grid points */
double meanFieldStrength(ref_ptr<VectorGrid> grid);

/** Evaluate the RMS of all grid points */
double rmsFieldStrength(ref_ptr<ScalarGrid> grid);
/** Evaluate the RMS of all grid points */
double rmsFieldStrength(ref_ptr<VectorGrid> grid);

/** Multiply all grid values by a given factor */
void scaleGrid(ref_ptr<ScalarGrid> grid, double a);
/** Multiply all grid values by a given factor */
void scaleGrid(ref_ptr<VectorGrid> grid, double a);

#ifdef CRPROPA_HAVE_FFTW3F
/**
 Create a random initialization of a turbulent field.
 @param lMin	Minimum wavelength of the turbulence
 @param lMax	Maximum wavelength of the turbulence
 @param alpha	Power law index of <B^2(k)> ~ k^alpha (alpha = -11/3 corresponds to a Kolmogorov spectrum)
 @param Brms	RMS field strength
 @param seed	Random seed
 */
void initTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin, double lMax, 
	   double alpha = -11./3., int seed = 0, bool helicity = false, double H = 0);
#endif // CRPROPA_HAVE_FFTW3F

/** Analytically calculate the correlation length of a turbulent field */
double turbulentCorrelationLength(double lMin, double lMax,
		double alpha = (-11./3.));

/** Fill vector grid from provided magnetic field */
void fromMagneticField(ref_ptr<VectorGrid> grid, ref_ptr<MagneticField> field);

/** Fill scalar grid from provided magnetic field */
void fromMagneticFieldStrength(ref_ptr<ScalarGrid> grid, ref_ptr<MagneticField> field);

/** Load a VectorGrid from a binary file with single precision */
void loadGrid(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Load a ScalarGrid from a binary file with single precision */
void loadGrid(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a VectorGrid to a binary file */
void dumpGrid(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a ScalarGrid to a binary file with single precision */
void dumpGrid(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** Load a VectorGrid grid from a plain text file */
void loadGridFromTxt(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Load a ScalarGrid from a plain text file */
void loadGridFromTxt(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a VectorGrid to a plain text file */
void dumpGridToTxt(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a ScalarGrid to a plain text file */
void dumpGridToTxt(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** @}*/
} // namespace crpropa

#endif // CRPROPA_GRIDTOOLS_H

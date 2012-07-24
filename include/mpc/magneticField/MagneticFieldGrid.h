#ifndef MPC_MAGNETICFIELDGRID_H_
#define MPC_MAGNETICFIELDGRID_H_

#include "mpc/magneticField/MagneticField.h"
#include "mpc/magneticField/PeriodicGrid.h"

namespace mpc {

/**
 @class MagneticFieldGrid
 @brief Magnetic field on a periodic, cartesian grid with trilinear interpolation.

 This class provides a container for a VectorFieldGrid (PeriodicGrid of type Vector3f) to serve as a MagneticField.
 */
class MagneticFieldGrid: public MagneticField {
	VectorFieldGrid *grid;

public:
	MagneticFieldGrid(VectorFieldGrid *grid);
	void setGrid(VectorFieldGrid *grid);
	VectorFieldGrid *getGrid();
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class MagneticFieldGrid
 @brief Modulated magnetic field on a periodic grid.

 The class
 */
class ModulatedMagneticFieldGrid: public MagneticField {
	VectorFieldGrid *grid;
	ScalarFieldGrid *modGrid;

public:
	ModulatedMagneticFieldGrid(VectorFieldGrid *grid, ScalarFieldGrid *modGrid);
	void setGrid(VectorFieldGrid *grid);
	void setModulationGrid(ScalarFieldGrid *modGrid);
	VectorFieldGrid *getGrid();
	ScalarFieldGrid *getModulationGrid();
	Vector3d getField(const Vector3d &position) const;
};

/** Calculate the mean field strength */
Vector3f meanFieldStrength(VectorFieldGrid *m);

/** Calculate the RMS field strength */
double rmsFieldStrength(VectorFieldGrid *m);

/** Multiply a magnetic field grid by a factor. */
void scale(VectorFieldGrid *m, double a);

#ifdef MPC_HAVE_FFTW3F
/**
 * Create a random initialization of a turbulent field.
 * @param lMin		Minimum wavelength of the turbulence
 * @param lMax		Maximum wavelength of the turbulence
 * @param spectralIndex	Power spectral index of the turbulence (-11/3 corresponds to a Kolmogorov spectrum)
 * @param Brms		RMS field strength
 * @param seed		Random seed
 */
void initTurbulence(VectorFieldGrid *m, double Brms, double lMin,
		double lMax, double spectralIndex = -11. / 3., int seed = 0);

/** Analytically calculate the correlation length of a turbulent field */
double turbulentCorrelationLength(double lMin, double lMax, double spectralIndex = -11. / 3.);
#endif // MPC_HAVE_FFTW3F

/**
 * Load a magnetic field grid from a binary file.
 * The field is stored single precision numbers with the field components in xyz order and the grid z-index changing the fastest.
 * @param c	conversion of the values in the file to SI
 */
void load(VectorFieldGrid *m, std::string filename, double c = 1);

/**
 * Dump a magnetic field grid to a binary file.
 * The field is stored single precision numbers with the field components in xyz order and the grid z-index changing the fastest.
 * @param c	conversion of the values in the file to SI
 */
void dump(VectorFieldGrid *m, std::string filename, double c = 1);

/**
 * Load a magnetic field grid from a plain text file.
 * The field is stored as one grid point per line from (0,0,0) to (nx, ny, nz) with the grid z-index changing the fastest.
 * Within one line the magnetic field values are stored in xyz order separated by a blank or tab.
 * Header lines must start with a #.
 * @param c	conversion of the values in the file to SI
 */
void loadTxt(VectorFieldGrid *m, std::string filename, double c = 1);

/**
 * Dump a magnetic field grid to a plain text file.
 * The field is stored as one grid point per line from (0,0,0) to (nx, ny, nz) with the grid z-index changing the fastest.
 * Within one line the magnetic field values are stored in xyz order separated by a blank or tab.
 * Header lines must start with a #.
 * @param c	conversion of SI to the values in the file
 */
void dumpTxt(VectorFieldGrid *m, std::string filename, double c = 1);

} // namespace mpc

#endif /* MPC_MAGNETICFIELDGRID_H_ */

#ifndef CRPROPA_GRIDTURBULENCE_H
#define CRPROPA_GRIDTURBULENCE_H

#ifdef CRPROPA_HAVE_FFTW3F

#include "crpropa/Grid.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"

#include "fftw3.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class GridTurbulence
 @brief Turbulent grid-based magnetic field with a general energy spectrum
 */
class GridTurbulence : public TurbulentField {
  protected:
	unsigned int seed;
	ref_ptr<Grid3f> gridPtr;

	void initGrid(const GridProperties &grid);
	void initTurbulence();

  public:
	/**
	 Create a random initialization of a turbulent field.
	 @param spectrum    TurbulenceSpectrum instance to define the spectrum of
	 turbulence
	 @param gridProp	GridProperties instance to define the underlying grid
	 @param seed	 Random seed
	 */
	GridTurbulence(const TurbulenceSpectrum &spectrum,
	               const GridProperties &gridProp, unsigned int seed = 0);

	Vector3d getField(const Vector3d &pos) const;

	/** Return a const reference to the grid */
	const ref_ptr<Grid3f> &getGrid() const;

	/* Helper functions for synthetic turbulent field models */
	// Check the grid properties before the FFT procedure
	static void checkGridRequirements(ref_ptr<Grid3f> grid, double lMin,
	                                  double lMax);
	// Execute inverse discrete FFT in-place for a 3D grid, from complex to real
	// space
	static void executeInverseFFTInplace(ref_ptr<Grid3f> grid,
	                                     fftwf_complex *Bkx, fftwf_complex *Bky,
	                                     fftwf_complex *Bkz);

	// Usefull checks for a grid field
	/** Evaluate the mean vector of all grid points */
	Vector3f getMeanFieldVector() const;
	/** Evaluate the mean of all grid points */
	double getMeanFieldStrength() const;
	/** Evaluate the RMS of all grid points */
	double getRmsFieldStrength() const;
	/** Evaluate the RMS of all grid points per axis */
	std::array<float, 3> getRmsFieldStrengthPerAxis() const;
	/** Evaluate generated power-spectrum */
	std::vector<std::pair<int, float>> getPowerSpectrum() const;
	/** Dump a Grid3f to a binary file */
	void dumpToFile(std::string filename) const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F

#endif // CRPROPA_GRIDTURBULENCE_H

#ifndef MPC_TURBULENTMAGNETICFIELD_H_
#define MPC_TURBULENTMAGNETICFIELD_H_

#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/Random.h"

namespace mpc {

/**
 @class TurbulentMagneticFieldGrid
 @brief Turbulent magnetic field on a cubic grid with trilinear interpolation.

 This class creates a random magnetic field with a turbulent spectrum.
 The field is isotropic and homogeneous with a zero mean magnetic field and and divergence.
 */
class TurbulentMagneticFieldGrid: public MagneticFieldGrid {
public:
	/**
	 * Constructor
	 * @param origin	Lower left front corner of the grid
	 * @param size 		Physical extension of the field grid
	 * @param samples	Number of grid samples per edge
	 */
	TurbulentMagneticFieldGrid(Vector3d origin, double size, size_t samples);

	/**
	 * Constructor, also performs initialization and normalization.
	 * @param origin	Lower left front corner of the grid
	 * @param size 		Physical extension of the field grid
	 * @param samples	Number of grid samples per edge
	 * @param lMin		Minimum wavelength of the turbulence
	 * @param lMax		Maximum wavelength of the turbulence
	 * @param spectralIndex	Power spectral index of the turbulence
	 * @param Brms		RMS field strength
	 */
	TurbulentMagneticFieldGrid(Vector3d origin, double size, size_t samples,
			double lMin, double lMax, double spectralIndex, double Brms);

	/**
	 * Define the properties of the turbulence.
	 * @param lMin		Minimum wavelength of the turbulence
	 * @param lMax		Maximum wavelength of the turbulence
	 * @param spectralIndex	Power spectral index turbulence
	 */
	void setTurbulenceProperties(double lMin, double lMax, double spectralIndex);

	/** Create a random initialization of the turbulent field. */
	void initialize();

	/** Create a random initialization from a given seed. */
	void initialize(int seed);

	double getPowerSpectralIndex() const;
	double getMinimumWavelength() const;
	double getMaximumWavelength() const;

	/** Return the analytical calculation of the correlation length. */
	double getCorrelationLength() const;

protected:
	double lMin; /**< Minimum wavelength of the turbulence */
	double lMax; /**< Maximum  wavelength of the turbulence */
	double spectralIndex; /**< Power spectral index of the turbulence */
	Random random; /**< Random number generator instance */
};

} // namespace mpc

#endif /* MPC_TURBULENTMAGNETICFIELD_H_ */

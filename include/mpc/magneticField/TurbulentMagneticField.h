#ifndef TURBULENTMAGNETICFIELD_H_
#define TURBULENTMAGNETICFIELD_H_

#include "mpc/magneticField/MagneticField.h"
#include "mpc/Random.h"

#include <vector>

namespace mpc {

/**
 @class TurbulentMagneticField
 @brief Turbulent magnetic field defined by a superposition of N random modes

 This class represents random magnetic field with a turbulent spectrum.
 The field is calculated at any point with double precision from a number of random modes.
 For reference see Giacinti 2011, DOI: 10.1016/j.astropartphys.2011.07.006
 Note that the normalization of the turbulent modes is not yet correct.
 Note that this method is slow: O(0.1 ms) on a 2.1 GHz Centrino
 */
class TurbulentMagneticField: public MagneticField {
public:
	TurbulentMagneticField() :
			nModes(0), spectralIndex(0), Brms(0), lMin(0), lMax(0) {
	}

	/** Constructor, also initializes the field. */
	TurbulentMagneticField(double Brms, double lMin, double lMax,
			double spectralIndex = -11. / 3., int nModes = 1000);

	/** Calculates the magnetic field at position from the dialed random turbulent modes */
	Vector3d getField(const Vector3d &position) const;

	/**
	 * Define the properties of the turbulence.
	 * @param Brms		RMS field strength
	 * @param lMin		Minimum wavelength of the turbulence
	 * @param lMax		Maximum wavelength of the turbulence
	 * @param spectralIndex	Power spectral index turbulence
	 * @param modes		Number of modes
	 */
	void setTurbulenceProperties(double Brms, double lMin, double lMax,
			double spectralIndex = -11. / 3., int nModes = 1000);

	/** Create a random initialization of the turbulent field. */
	void initialize();

	/** Create a random initialization from a given seed. */
	void initialize(int seed);

	double getPowerSpectralIndex() const;
	double getMinimumWavelength() const;
	double getMaximumWavelength() const;
	double getRMSFieldStrength() const;

	/** Return the analytical calculation of the correlation length. */
	double getCorrelationLength() const;

private:
	struct Mode {
		double amplitude;
		double phase;
		Vector3d e1;
		Vector3d e2;
		Vector3d k;
	};
	std::vector<Mode> modes; /**< List of random turbulent modes */
	int nModes; /**< Number of modes */
	double spectralIndex; /**< Power spectral index of the turbulence */
	double Brms; /**< RMS Field strength */
	double lMin; /**< Minimum wavelength */
	double lMax; /**< Maximum wavelength */
	Random random; /**< Random number generator instance */
};

} // mpc

#endif /* TURBULENTMAGNETICFIELD_H_ */

#ifndef MPC_SPHTURBULENTMAGNETICFIELDGRID_H_
#define MPC_SPHTURBULENTMAGNETICFIELDGRID_H_

#ifdef MPC_HAVE_GADGET

#include "mpc/magneticField/TurbulentMagneticFieldGrid.h"
#include "mpc/magneticField/SPHMagneticField.h"

#include <memory>

namespace mpc {

/**
 @class SPHTurbulentMagneticFieldGrid
 @brief Random turbulent magnetic field on a cubic grid modulated with an SPH density field.
 */
class SPHTurbulentMagneticFieldGrid: public TurbulentMagneticFieldGrid {
public:
	/** Constructor. Reimplementation of TurbulentMagneticField. */
	SPHTurbulentMagneticFieldGrid(Vector3d origin, double size, size_t samples) :
			TurbulentMagneticFieldGrid(origin, size, samples) {
	}

	/** Constructor. Reimplementation of TurbulentMagneticField. */
	SPHTurbulentMagneticFieldGrid(Vector3d origin, double size, size_t samples,
			double lMin, double lMax, double spectralIndex, double Brms) :
			TurbulentMagneticFieldGrid(origin, size, samples, lMin, lMax,
					spectralIndex, Brms) {
	}

	/**
	 * Modulate the magnetic field with the baryon density from a gadget SPH field.
	 * @param filename	Path to the SPH database file
	 * @param origin	Origin of the SPH field
	 * @param size		Size of the SPH field
	 * @param bins		Number of bins of the SPH field
	 * @param exp		Exponent for the density modulation
	 * @param norm		Normalization of the field strength
	 */
	void setModulation(std::string filename, Vector3d origin, double size,
			size_t samples, double exp, double norm);

	/** Get the density modulated magnetic field vector */
	Vector3d getField(const Vector3d &position) const;

	/** Get the density in [10^10 M_sol * h^2 / kpc^3] */
	double getRho(const Vector3d &position) const;

	double exponent; /**< Exponent of modulation with the density field */
	double normalization; /**< Normalization factor */
	std::auto_ptr<SPHMagneticField> sphField; /**< SPH (density) field */
};

} // namespace mpc

#endif // MPC_HAVE_GADGET
#endif // MPC_SPHTURBULENTMAGNETICFIELDGRID_H_

#ifndef MPC_SPHTURBULENTMAGNETICFIELD_H_
#define MPC_SPHTURBULENTMAGNETICFIELD_H_

#include "mpc/magneticField/TurbulentMagneticField.h"

#ifdef MPC_HAVE_GADGET

namespace mpc {

/**
 @class SPHTurbulentMagneticField
 @brief Random turbulent magnetic field on a cubic grid modulated with large scale structure.
 */
class SPHTurbulentMagneticField: public TurbulentMagneticField {
public:
	/** Constructor. Reimplementation of TurbulentMagneticField. */
	SPHTurbulentMagneticField(Vector3d origin, double size, size_t samples) :
			TurbulentMagneticField(origin, size, samples) {
	}

	/** Constructor. Reimplementation of TurbulentMagneticField. */
	SPHTurbulentMagneticField(Vector3d origin, double size, size_t samples,
			double lMin, double lMax, double spectralIndex, double Brms) :
			TurbulentMagneticField(origin, size, samples, lMin, lMax,
					spectralIndex, Brms) {
	}

	/**
	 * Modulate the magnetic field with the baryon density from a gadget SPH field.
	 * The field should be (re)normalized afterwards.
	 * @param filename	Path to the SPH database file
	 * @param exp		Exponent for the modulation
	 */
	void modulate(std::string filename, double exp);
};

} // namespace mpc

#endif // MPC_HAVE_GADGET
#endif // MPC_SPHTURBULENTMAGNETICFIELD_H_

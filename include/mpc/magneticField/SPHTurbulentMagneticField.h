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
	SPHTurbulentMagneticField(Vector3d origin, size_t samples, double spacing) :
			TurbulentMagneticField(origin, samples, spacing) {
	}
	void modulate(std::string filename);
};

} // namespace mpc

#endif // MPC_HAVE_GADGET
#endif // MPC_SPHTURBULENTMAGNETICFIELD_H_

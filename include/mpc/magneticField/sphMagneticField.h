#ifndef SPHMAGNETICFIELD_H_
#define SPHMAGNETICFIELD_H_

#include "mpc/magneticField/magneticFieldGrid.h"

#include "gadget/MagneticField.h"
#include "gadget/SmoothParticle.h"
#include "gadget/Vector3.h"

#include <vector>
#include <memory>

namespace mpc {

/**
 @class SPHMagneticField
 @brief Wrapper for gadget::DirectMagneticField
 */
class SPHMagneticField: public MagneticField {
	size_t gridSize;
public:
	SPHMagneticField(Vector3 origin, double size, size_t gridSize);
	Vector3 getField(const Vector3 &position) const;
	void updateSimulationVolume(const Vector3 &origin, double size);

	std::auto_ptr<gadget::DirectMagneticField> gadgetField;
};

/**
 @class SPHMagneticFieldGrid
 @brief Wrapper for gadget::SampledMagneticField
 */
class SPHMagneticFieldGrid: public MagneticField {
public:
	SPHMagneticFieldGrid(Vector3 origin, double size, size_t n);
	Vector3 getField(const Vector3 &position) const;
	void updateSimulationVolume(const Vector3 &origin, double size);

	std::auto_ptr<gadget::SampledMagneticField> gadgetField;
};

} // namespace mpc

#endif /* SPHMAGNETICFIELD_H_ */

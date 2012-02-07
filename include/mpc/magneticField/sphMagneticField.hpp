#ifndef SPHMAGNETICFIELD_HPP_
#define SPHMAGNETICFIELD_HPP_

#include "mpc/magneticField/magneticFieldGrid.hpp"
#include "gadget/MagneticField.hpp"
#include "gadget/SmoothParticle.hpp"
#include "gadget/Vector3.hpp"
#include <vector>
#include <memory>

namespace mpc {

/**
 @class SPHMagneticField
 @brief Wrapper for gadget::DirectMagneticField
 */
class SPHMagneticField: public MagneticField {
public:
	SPHMagneticField(Vector3 origin, double size, size_t n);
	Vector3 getField(const Vector3 &position) const;

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

	std::auto_ptr<gadget::SampledMagneticField> gadgetField;
};

} // namespace mpc

#endif /* SPHMAGNETICFIELD_HPP_ */

#ifndef MPC_SPHMAGNETICFIELD_H_
#define MPC_SPHMAGNETICFIELD_H_

#ifdef MPC_HAVE_GADGET

#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/Units.h"

#include "gadget/Database.h"
#include "gadget/MagneticField.h"

#include <vector>
#include <memory>

namespace gadget {
	class DirectMagneticField;
	class SampledMagneticField;
}

namespace mpc {

/**
 @class SPHMagneticField
 @brief Wrapper for gadget::DirectMagneticField
 */
class SPHMagneticField: public MagneticField {
	gadget::DirectMagneticField field;
	gadget::FileDatabase database;
public:
	SPHMagneticField(Vector3d origin, double size, size_t samples,
			std::string filename);
	SPHMagneticField(size_t samples, std::string filename);

	Vector3d getField(const Vector3d &position) const;
	double getRho(const Vector3d &position) const; /*< Return the density in [kg / m^3] */
};

/**
 @class SPHMagneticFieldGrid
 @brief Wrapper for gadget::SampledMagneticField
 */
class SPHMagneticFieldGrid: public MagneticField {
	gadget::SampledMagneticField field;
	gadget::FileDatabase database;
	size_t samples;
	std::string cachePrefix;
	bool cacheEnabled;
public:
	SPHMagneticFieldGrid(Vector3d origin, double size, size_t gridSize,
			std::string filename);
	SPHMagneticFieldGrid(size_t gridSize, std::string filename);
	void init(const Vector3d &origin, double size);
	Vector3d getField(const Vector3d &position) const;
	void setCachePrefix(std::string prefix);
	void setCacheEnabled(bool enabled);
};

/**
 @class GadgetMagneticField
 @brief Wrapper for gadget::MagneticField
 */
class GadgetMagneticField: public MagneticField {
	gadget::ref_ptr<gadget::MagneticField> field;
public:
	GadgetMagneticField(gadget::ref_ptr<gadget::MagneticField> field) : field(field) {

	}
	Vector3d getField(const Vector3d &position) const {
		gadget::Vector3f b, r = gadget::Vector3f(position.x, position.y, position.z);
		bool isGood = field->getField(r / kpc, b);
		if (!isGood)
			std::cerr << "mpc::SPHMagneticField invalid position : " << position
					<< std::endl;
		return Vector3d(b.x, b.y, b.z) * gauss;
	}
};

} // namespace mpc

#endif // MPC_HAVE_GADGET
#endif // MPC_SPHMAGNETICFIELD_H_

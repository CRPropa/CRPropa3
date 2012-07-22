#ifndef MPC_SPHMAGNETICFIELD_H_
#define MPC_SPHMAGNETICFIELD_H_

#ifdef MPC_HAVE_GADGET

#include "mpc/magneticField/MagneticFieldGrid.h"

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
	Vector3d getField(const Vector3d &position) const;
	void setCachePrefix(std::string prefix);
	void setCacheEnabled(bool enabled);
};

} // namespace mpc

#endif // MPC_HAVE_GADGET
#endif // MPC_SPHMAGNETICFIELD_H_

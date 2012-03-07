#ifndef SPHMAGNETICFIELD_H_
#define SPHMAGNETICFIELD_H_

#include "mpc/magneticField/magneticFieldGrid.h"

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
	SPHMagneticField(Vector3 origin, double size, size_t gridSize,
			const std::string filename);
	Vector3 getField(const Vector3 &position) const;
	void updateSimulationVolume(const Vector3 &origin, double size);

};

/**
 @class SPHMagneticFieldGrid
 @brief Wrapper for gadget::SampledMagneticField
 */
class SPHMagneticFieldGrid: public MagneticField {
	gadget::SampledMagneticField field;
	gadget::FileDatabase database;
public:
	SPHMagneticFieldGrid(Vector3 origin, double size, size_t samples,
			const std::string filename);
	Vector3 getField(const Vector3 &position) const;
	void updateSimulationVolume(const Vector3 &origin, double size);

};

} // namespace mpc

#endif /* SPHMAGNETICFIELD_H_ */

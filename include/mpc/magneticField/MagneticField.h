#ifndef MPC_MAGNETICFIELD_H_
#define MPC_MAGNETICFIELD_H_

#include "mpc/Vector3.h"
#include "mpc/Referenced.h"

namespace mpc {

/**
 @class MagneticField
 @brief Abstract base class for magnetic fields.
 */
class MagneticField: public Referenced {
public:
	virtual ~MagneticField() {
	}
	virtual Vector3d getField(const Vector3d &position) const = 0;
};

/**
 @class PeriodicMagneticField
 @brief Magnetic field decorator to implement periodic magnetic fields.
 */
class PeriodicMagneticField: public MagneticField {
	ref_ptr<MagneticField> field;
	Vector3d origin, extends;
	bool reflective;
public:
	PeriodicMagneticField(ref_ptr<MagneticField> field,
			const Vector3d &extends);
	PeriodicMagneticField(ref_ptr<MagneticField> field, const Vector3d &extends,
			const Vector3d &origin, bool reflective);
	Vector3d &getOrigin();
	void setOrigin(const Vector3d &origin);
	Vector3d &getExtends();
	void setExtends(const Vector3d &origin);
	bool isReflective();
	void setReflective(bool reflective);
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class MagneticFieldList
 @brief List of magnetic fields

 The field at a given position is the sum of all fields evaluated at that position.
 */
class MagneticFieldList: public MagneticField {
	std::vector<ref_ptr<MagneticField> > fields;
public:
	void addField(ref_ptr<MagneticField> field);
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class UniformMagneticField
 @brief Magnetic field with one B-field vector.
 */
class UniformMagneticField: public MagneticField {
private:
	Vector3d value;

public:
	UniformMagneticField(const Vector3d &value) :
			value(value) {
	}
	Vector3d getField(const Vector3d &position) const {
		return value;
	}
};

} // namespace mpc

#endif /* MPC_MAGNETICFIELD_H_ */

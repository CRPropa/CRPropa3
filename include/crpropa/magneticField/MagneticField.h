#ifndef CRPROPA_MAGNETICFIELD_H
#define CRPROPA_MAGNETICFIELD_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {

/**
 @class MagneticField
 @brief Abstract base class for magnetic fields.
 */
class MagneticField: public Referenced {
public:
	virtual ~MagneticField() {
	}
	virtual Vector3d getField(const Vector3d &position) const {};
	virtual Vector3d getField(const Vector3d &position, double z) const {
		return getField(position);
	};
};

/**
 @class PeriodicMagneticField
 @brief Magnetic field decorator implementing periodic fields.
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
 @brief Magnetic field decorator implementing a superposition of fields.
 */
class MagneticFieldList: public MagneticField {
	std::vector<ref_ptr<MagneticField> > fields;
public:
	void addField(ref_ptr<MagneticField> field);
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class MagneticFieldEvolution
 @brief Magnetic field decorator implementing an evolution of type (1+z)^m.
 */
class MagneticFieldEvolution: public MagneticField {
	ref_ptr<MagneticField> field;
	double m;
public:
	MagneticFieldEvolution(ref_ptr<MagneticField> field, double m);
	Vector3d getField(const Vector3d &position, double z = 0) const;
};

/**
 @class UniformMagneticField
 @brief Magnetic field with one B-field vector.
 */
class UniformMagneticField: public MagneticField {
	Vector3d value;
public:
	UniformMagneticField(const Vector3d &value) :
			value(value) {
	}
	Vector3d getField(const Vector3d &position) const {
		return value;
	}
};

/**
 @class MagneticDipoleField
 @brief Magnetic dipole field defined by the magnetic moment and the 'core' radius.
 */
class MagneticDipoleField: public MagneticField {
	Vector3d origin;
	Vector3d moment;
	double radius;
public:
	MagneticDipoleField(const Vector3d &pos, const Vector3d &moment, const double radius) :
			origin(origin), moment(moment), radius(radius) {
	}
	Vector3d getField(const Vector3d &position) const {
		Vector3d r = (position - origin);
		Vector3d unit_r = r.getUnitVector();
		
		if (r.getR() == 0) // skip singularity
			return Vector3d(0, 0, 0);
		return unit_r * moment.dot(unit_r) / pow((r.getR()/radius), 3) * mu0 / (4*M_PI);
	}
};

} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELD_H

#ifndef CRPROPA_MAGNETICFIELD_H
#define CRPROPA_MAGNETICFIELD_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

#ifdef CRPROPA_HAVE_MUPARSER
#include "muParser.h"
#endif

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class MagneticField
 @brief Abstract base class for magnetic fields.
 */
class MagneticField: public Referenced {
public:
	virtual ~MagneticField() {
	}
	virtual Vector3d getField(const Vector3d &position) const {
		return Vector3d(0,0,0);
	};
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
	MagneticDipoleField(const Vector3d &origin, const Vector3d &moment, const double radius) :
			origin(origin), moment(moment), radius(radius) {
	}
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class SingleModeHelicalMagneticField
 @brief helical magnetic field with a single polarization mode defined by
 - the wavelength (wavelength) of the polarization mode,
 - its handedness (handedness),
 - the magnetic field amplitude at the origin (amplitudeOrigin),
 - a unit vector pointing in the direction of the magnetic field at the origin (unitVectorOrigin), and
 - a second unit vector (unitVector2) which, together with unitVectorOrigin, spans the polarization plane.
**/
class SingleModeHelicalMagneticField: public MagneticField {
	Vector3d origin;
	Vector3d unitVectorOrigin;
	Vector3d unitVector2;
	double amplitudeOrigin;
	double wavelength;
	double handedness;
public:
	SingleModeHelicalMagneticField(const Vector3d &origin, const Vector3d &unitVectorOrigin, const Vector3d &unitVector2, const double amplitudeOrigin, const double wavelength, const double handedness) :
			origin(origin),  unitVectorOrigin(unitVectorOrigin),  unitVector2(unitVector2),  amplitudeOrigin(amplitudeOrigin), wavelength(wavelength), handedness(handedness) {
	}
	Vector3d getField(const Vector3d &position) const;
};

#ifdef CRPROPA_HAVE_MUPARSER
/**
 @class RenormalizeMagneticField
 @brief Renormalize strength of a given field by expression in which B is the strength variable.
 */
class RenormalizeMagneticField: public MagneticField {
	ref_ptr<MagneticField> field;
	std::string expression;
	mu::Parser *p;
	double Bmag;
public:
	RenormalizeMagneticField(ref_ptr<MagneticField> field, std::string expression);
	~RenormalizeMagneticField() { delete p;	}
	Vector3d getField(const Vector3d &position);
};
#endif

/** @} */
} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELD_H

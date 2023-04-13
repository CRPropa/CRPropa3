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

 The periodic cube is defined by its origin (Vector3d) and an
 extends parameter (Vector3d). All points x=(x_1, x_2, x_3) 
 that are described by x_i = origin_i + epsilon * extend_i, 
 with epsilon = 0...1 are within the base cube. Magnetic field
 strengths for all positions outside of this cube are calculated 
 based on the values in the base cube. 
 This can be done periodically or reflectively.
*/

class PeriodicMagneticField: public MagneticField {
	ref_ptr<MagneticField> field;
	Vector3d origin, extends;
	bool reflective;
public:
	/**
	 * Constructor
	 * @param field magnetic field reference pointer
	 * @param extends length, width, and height of the base cube 
	*/
	PeriodicMagneticField(ref_ptr<MagneticField> field,
			const Vector3d &extends);
	/**
	 * Constructor
	 * @param field magnetic field reference pointer
	 * @param extends length, width, and height of the base cube
	 * @param origin defines the reference position 
	 * @param reflective for periodic or reflective behavior  
	*/
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
	/**
	 * Constructor
	 * @param field magnetic field reference pointer
	 * @param m cosmic evolution parameter 
	*/
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
	/**
	 * Constructor
	 * @param value magnetic field strength
	*/
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
	/**
	 * Constructor
	 * @param origin 	singularity of the dipole field
	 * @param moment 	magnetic moment of the dipole field
	 * @param radius 	inside a radius around the origin the 
	 * 					magnetic field is constant: moment * 2 * mu0 / 3  
	*/
	MagneticDipoleField(const Vector3d &origin, const Vector3d &moment, const double radius) :
			origin(origin), moment(moment), radius(radius) {
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
	/**
	 * Constructor
	 * @param field 		magnetic field reference pointer
	 * @param expression 	muParser expression used to renormalize the field, 
	 * 						e.g., "gauss". 
	*/
	RenormalizeMagneticField(ref_ptr<MagneticField> field, std::string expression);
	~RenormalizeMagneticField() { delete p;	}
	Vector3d getField(const Vector3d &position);
};
#endif

/** @} */
} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELD_H

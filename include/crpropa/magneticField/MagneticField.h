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
 @class MagneticBottle
 @brief Magnetic bottle along the z-axis. The field has the constant value Bz_strength for z in [-width, width]. 
		Out of this range, the field lines run together and the B-field strength raises linearly.
		This B-field was inspired by the magnetic bottle in Computer Algebra Recipes by Enns and McGuire
 **/
class MagneticBottle: public MagneticField {
	
	double Bz_max;
	double Bz_min;
	double L;
	
public:
	MagneticBottle(double Bz_max, double Bz_min, double lenght) : 
			 Bz_max(Bz_max), Bz_min(Bz_min), L(lenght) {}
	Vector3d getField(const Vector3d &position) const;	
	
};

class GyroField: public MagneticField {

	const double B_max;
	const double B_min;
	const double radius;
	
public:
	GyroField(const double B_max, const double B_min, const double radius) : 
			 B_max(B_max), B_min(B_min), radius(radius) {}
	Vector3d getField(const Vector3d &position) const;	
	
};

class LongConductorField: public MagneticField {

	const double B_radius;
	const double radius;
	
public:
	LongConductorField(const double B_radius, const double radius) : 
			 B_radius(B_radius), radius(radius) {}
	Vector3d getField(const Vector3d &position) const;	
	
};

class CircleField: public MagneticField {

	const double B;
	
public:
	CircleField(const double B) : 
			 B(B) {}
	Vector3d getField(const Vector3d &position) const;	
	
};

class HongQinField: public MagneticField {

	const double B;
	const double delta;
	
public:
	HongQinField(const double B, const double delta) : 
			 B(B), delta(delta) {}
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

} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELD_H

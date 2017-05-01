#ifndef CRPROPA_ADVECTIONFIELD_H
#define CRPROPA_ADVECTIONFIELD_H

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {

/**
 @class AdvectionField
 @brief Abstract base class for advection fields. These are used to model 
	the deterministic part of the Fokker-Planck equation.
 */
class AdvectionField: public Referenced {
public:
	virtual ~AdvectionField() {
	}
	virtual Vector3d getField(const Vector3d &position) const {};
};


/**
 @class AdvectionFieldList
 @brief Advection field decorator implementing a superposition of fields.
 */
class AdvectionFieldList: public AdvectionField {
	std::vector<ref_ptr<AdvectionField> > fields;
public:
	void addField(ref_ptr<AdvectionField> field);
	Vector3d getField(const Vector3d &position) const;
};


/**
 @class UniformAdvectionField
 @brief Advection field with one B-field vector.
 */
class UniformAdvectionField: public AdvectionField {
	Vector3d value;
public:
	UniformAdvectionField(const Vector3d &value);
	Vector3d getField(const Vector3d &position) const;
};

} // namespace crpropa

#endif // CRPROPA_ADVECTIONFIELD_H

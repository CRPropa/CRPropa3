#ifndef CRPROPA_ADVECTIONFIELD_H
#define CRPROPA_ADVECTIONFIELD_H

#pragma once
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"

namespace crpropa {

/**
 @class AdvectionField
 @brief Abstract base class for advection fields. These are used to model 
	the deterministic part of the Fokker-Planck equation. The getDivergence()
	method is used to model the adibatic cooling/heating.
 */
class AdvectionField: public Referenced {
public:
	virtual ~AdvectionField() {
	}
	virtual Vector3d getField(const Vector3d &position) const {};
	virtual double getDivergence(const Vector3d &position) const {};
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
	double getDivergence(const Vector3d &position) const;
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
	double getDivergence(const Vector3d &position) const;
};


/**
 @class SphericalAdvectionField
 @brief Spherical advection with a exponentially increasing and
	exponentially constant velocity.
*/

class SphericalAdvectionField: public AdvectionField {
	Vector3d origin; //origin of the advection sphere
	double radius; //radius of the advection sphere
	double vMax; // maximum wind velocity
	double tau; // transition distance
	double alpha; //tuning parameter
public:
	SphericalAdvectionField(Vector3d origin, double radius, double vMax, double tau, double alpha);
	Vector3d getField(const Vector3d &position) const;
	double getDivergence(const Vector3d &position) const;

	double getV(const double &r) const;

	void setOrigin(Vector3d origin);
	void setRadius(double radius);
	void setVMax(double vMax);
	void setTau(double tau);
	void setAlpha(double alpha);

	Vector3d getOrigin() const;
	double getRadius() const;
	double getVMax() const;
	double getTau() const;
	double getAlpha() const;
	
	std::string getDescription() const;
};


} // namespace crpropa

#endif // CRPROPA_ADVECTIONFIELD_H

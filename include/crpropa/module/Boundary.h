#ifndef CRPROPA_BOUNDARY_H
#define CRPROPA_BOUNDARY_H

#include "crpropa/Module.h"

namespace crpropa {
/**
 * \addtogroup Condition
 * @{
 */

/**
 @class PeriodicBox
 @brief Rectangular box with periodic boundaries.

 If a particle passes on of the sides it is placed at the opposite side and its initial (source) position changed accordingly.
 This implements periodic boundaries, that keep the particle inside the box and instead move the source away periodically.
 Particles can overshoot (be outside of the box during the step) since the step size is not limited by this module.
 */
class PeriodicBox: public Module {
private:
	Vector3d origin;
	Vector3d size;

public:
	PeriodicBox();
	PeriodicBox(Vector3d origin, Vector3d size);
	void process(Candidate *candidate) const;
	void setOrigin(Vector3d origin);
	void setSize(Vector3d size);
	std::string getDescription() const;
};

/**
 @class ReflectiveBox
 @brief Rectangular box with reflective boundaries.

 If a particle passes on of the sides it is reflected back inside (position and velocity) and its initial position changed as if the particle had come from that side.
 This implements periodic boundaries, that keep the particle inside the box and instead move the source away reflectively.
 Particles can overshoot (be outside of the box during the step) since the step size is not limited by this module.
 */
class ReflectiveBox: public Module {
private:
	Vector3d origin;
	Vector3d size;

public:
	ReflectiveBox();
	ReflectiveBox(Vector3d origin, Vector3d size);
	void process(Candidate *candidate) const;
	void setOrigin(Vector3d origin);
	void setSize(Vector3d size);
	std::string getDescription() const;
};

/**
 @class CubicBoundary
 @brief Flags a particle when exiting the cube.

 The particle is made inactive and flagged as "Rejected".
 By default the module prevents overshooting the boundary by more than a margin of 0.1 kpc.
 This corresponds to the default minimum step size of the propagation modules (PropagationCK and SimplePropagation).
 */
class CubicBoundary: public AbstractCondition {
private:
	Vector3d origin;
	double size;
	double margin;
	bool limitStep;

public:
	CubicBoundary();
	CubicBoundary(Vector3d origin, double size);
	void process(Candidate *candidate) const;
	void setOrigin(Vector3d origin);
	void setSize(double size);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	std::string getDescription() const;
};

/**
 @class SphericalBoundary
 @brief Flag a particle when leaving the sphere.

 The particle is made inactive and flagged as "Rejected".
 By default the module prevents overshooting the boundary by more than a margin of 0.1 kpc.
 This corresponds to the default minimum step size of the propagation modules (PropagationCK and SimplePropagation).
 */
class SphericalBoundary: public AbstractCondition {
private:
	Vector3d center;
	double radius;
	double margin;
	bool limitStep;

public:
	SphericalBoundary();
	SphericalBoundary(Vector3d center, double radius);
	void process(Candidate *candidate) const;
	void setCenter(Vector3d center);
	void setRadius(double size);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	std::string getDescription() const;
};

/**
 @class EllipsoidalBoundary
 @brief Flags a particle when leaving the ellipsoid.

 This module flags particles when outside of the ellipsoid, defined by two focal points and a major axis (length).
 The particle is made inactive and flagged as "Rejected".
 By default the module prevents overshooting the boundary by more than a margin of 0.1 kpc.
 This corresponds to the default minimum step size of the propagation modules (PropagationCK and SimplePropagation).
 */
class EllipsoidalBoundary: public AbstractCondition {
private:
	Vector3d focalPoint1;
	Vector3d focalPoint2;
	double majorAxis;
	double margin;
	bool limitStep;

public:
	EllipsoidalBoundary();
	EllipsoidalBoundary(Vector3d focalPoint1, Vector3d focalPoint2,
			double majorAxis);
	void process(Candidate *candidate) const;
	void setFocalPoints(Vector3d focalPoint1, Vector3d focalPoint2);
	void setMajorAxis(double size);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	std::string getDescription() const;
};


/**
 @class CylindricalBoundary
 @brief Flags a particle when leaving the cylinder.
 This module flags particles when outside of the cylinder, defined by a radius and a height.
 The particle is made inactive and by default is flagged "OutOfBounds".
 Optionally the module can ensure the candidate does not overshoot the boundary by more than a set margin.
 */
class CylindricalBoundary: public AbstractCondition {
private:
	Vector3d origin;
	double height;
	double radius;
double margin;
	bool limitStep;

public:
	CylindricalBoundary();
	CylindricalBoundary(Vector3d origin, double height,
			double radius);
	void process(Candidate *candidate) const;
	void setOrigin(Vector3d origin);
	void setHeight(double height);
	void setRadius(double radius);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_BOUNDARY_H

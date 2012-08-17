#ifndef MPC_BOUNDARY_H_
#define MPC_BOUNDARY_H_

#include "mpc/Module.h"

namespace mpc {

/**
 @class PeriodicBox
 @brief Rectangular box with periodic boundaries.

 If a particle passes on of the sides it is placed at the opposite side and its initial (source) position changed accordingly.
 This realizes periodic boundaries, where the particle is kept inside the box and the source is moved away from the box.
 Particles can overshoot (be outside of the box during the step) since the step size is not limited by this module.
 */
class PeriodicBox: public Module {
private:
	Vector3d origin;
	Vector3d size;
	void updateDescription();

public:
	PeriodicBox(Vector3d origin, Vector3d size);
	void process(Candidate *candidate) const;
};

/**
 @class ReflectiveBox
 @brief Rectangular box with reflective boundaries.

 If a particle passes on of the sides it is reflected back inside (position and velocity) and its initial position changed as if the particle had come from that side.
 This realizes reflective boundaries, where the particle is kept inside the box and the source is moved away from the box.
 Particles can overshoot (be outside of the box during the step) since the step size is not limited by this module.
 */
class ReflectiveBox: public Module {
private:
	Vector3d origin;
	Vector3d size;
	void updateDescription();

public:
	ReflectiveBox(Vector3d origin, Vector3d size);
	void process(Candidate *candidate) const;
};

/**
 @class CubicBoundary
 @brief Flags a particle when exiting the cube.

 This module flags particles when outside of the cube, defined by a lower corner and edge length.
 The particle is made inactive and by default is flagged "OutOfBounds".
 Optionally the module can ensure the candidate does not overshoot the boundary by more than a set margin.
 */
class CubicBoundary: public Module {
private:
	Vector3d origin;
	double size;
	double margin;
	std::string flag;
	std::string flagValue;
	bool limitStep;
	void updateDescription();

public:
	CubicBoundary(Vector3d origin, double size,
			std::string flag = "OutOfBounds", std::string flagValue = "");
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

/**
 @class SphericalBoundary
 @brief Flag a particle when leaving the sphere.

 This module flags particles when outside of the sphere, defined by a center and radius.
 The particle is made inactive and by default is flagged "OutOfBounds".
 Optionally the module can ensure the candidate does not overshoot the boundary by more than a set margin.
 */
class SphericalBoundary: public Module {
private:
	Vector3d center;
	double radius;
	double margin;
	std::string flag;
	std::string flagValue;
	bool limitStep;
	void updateDescription();

public:
	SphericalBoundary(Vector3d center, double radius, std::string flag =
			"OutOfBounds", std::string flagValue = "");
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

/**
 @class EllipsoidalBoundary
 @brief Flags a particle when leaving the ellipsoid.

 This module flags particles when outside of the ellipsoid, defined by two focal points and a major axis (length).
 The particle is made inactive and by default is flagged "OutOfBounds".
 Optionally the module can ensure the candidate does not overshoot the boundary by more than a set margin.
 */
class EllipsoidalBoundary: public Module {
private:
	Vector3d focalPoint1;
	Vector3d focalPoint2;
	double majorAxis;
	double margin;
	std::string flag;
	std::string flagValue;
	bool limitStep;
	void updateDescription();

public:
	EllipsoidalBoundary(Vector3d focalPoint1, Vector3d focalPoint2,
			double majorAxis, std::string flag = "OutOfBounds",
			std::string flagValue = "");
	void setLimitStep(bool limitStep, double margin);
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* MPC_BOUNDARY_H_ */

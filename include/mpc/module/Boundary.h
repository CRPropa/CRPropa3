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
 This realizes reflective boundaries, where the particle is kept inside the box and the source is moved away from the box.
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

public:
	CubicBoundary();
	CubicBoundary(Vector3d origin, double size);
	void process(Candidate *candidate) const;
	void setOrigin(Vector3d origin);
	void setSize(double size);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	void setFlag(std::string flag, std::string flagValue);
	std::string getDescription() const;
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

public:
	SphericalBoundary();
	SphericalBoundary(Vector3d center, double radius);
	void process(Candidate *candidate) const;
	void setCenter(Vector3d center);
	void setRadius(double size);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	void setFlag(std::string flag, std::string flagValue);
	std::string getDescription() const;
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

public:
	EllipsoidalBoundary();
	EllipsoidalBoundary(Vector3d focalPoint1, Vector3d focalPoint2,
			double majorAxis);
	void process(Candidate *candidate) const;
	void setFocalPoints(Vector3d focalPoint1, Vector3d focalPoint2);
	void setMajorAxis(double size);
	void setMargin(double margin);
	void setLimitStep(bool limitStep);
	void setFlag(std::string flag, std::string flagValue);
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_BOUNDARY_H_ */

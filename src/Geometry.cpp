#include <limits>
#include <cmath>
#include "kiss/logger.h"
#include "crpropa/Geometry.h"

namespace crpropa
{
Plane::Plane(const Vector3d& _x0, const Vector3d& _n) : x0(_x0), n(_n) {};

Plane::Plane(const Vector3d& _x0, const Vector3d& v1,const Vector3d& v2) : x0(_x0), n(0,0,0)
{
	n = v1.cross(v2);
	n /= n.getR();
};

double Plane::distance(const Vector3d &x) const
{
	Vector3d dX = x - x0;
	return n.dot(dX);
};


Sphere::Sphere(const Vector3d& _center, double _radius) : center(_center), radius(_radius) {};
double Sphere::distance(const Vector3d &point) const
{
	Vector3d dR = point - center;
	return dR.getR() - radius;
}


ParaxialBox::ParaxialBox(const Vector3d& _corner, const Vector3d& _size) : corner(_corner), size(_size) {};
double ParaxialBox::distance(const Vector3d &point) const
{
	Vector3d X = point - corner - size/2.;

	if ((fabs(X.x) <= size.x/2.) and (fabs(X.y) <= size.y/2.) and (fabs(X.z) <= size.z/2.))
	{ // inside the cube
		Vector3d Xp = size/2. - X.abs();
		double d = Xp.min();
		return -1. * d;
	}

	double a = std::max(0., fabs(X.x) - size.x/2.);
	double b = std::max(0., fabs(X.y) - size.y/2.);
	double c = std::max(0., fabs(X.z) - size.z/2.);

	return sqrt(a*a + b*b +c*c);
}


} // namespace

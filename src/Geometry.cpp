#include <limits>
#include <cmath>
#include "kiss/logger.h"
#include "crpropa/Geometry.h"

#include <iostream>
namespace crpropa
{
// Plane ------------------------------------------------------------------
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

std::string Plane::getDescription() const
{
	std::stringstream ss;
	ss << "Plane: " << std::endl
		 << "   x0: " << x0 << std::endl
		 << "    n: " << n << std::endl;
	return ss.str();
};

Vector3d Plane::normal(const Vector3d& point) const
{
  return n;
}


// Sphere ------------------------------------------------------------------
Sphere::Sphere(const Vector3d& _center, double _radius) : center(_center), radius(_radius) {};

double Sphere::distance(const Vector3d &point) const
{
	Vector3d dR = point - center;
	return dR.getR() - radius;
}

Vector3d Sphere::normal(const Vector3d& point) const
{
  Vector3d d = point-center;
  return d.getUnitVector();
}

std::string Sphere::getDescription() const
{
	std::stringstream ss;
	ss << "Sphere: " << std::endl
		 << "   Center: " << center << std::endl
		 << "   Radius: " << radius << std::endl;
	return ss.str();
};


// ParaxialBox -------------------------------------------------------------
ParaxialBox::ParaxialBox(const Vector3d& _corner, const Vector3d& _size) : corner(_corner), size(_size) {};
double ParaxialBox::distance(const Vector3d &point) const
{
	Vector3d X = point - corner - size/2.;

	if ((fabs(X.x) <= size.x/2.) and (fabs(X.y) <= size.y/2.) and (fabs(X.z) <= size.z/2.))
	{ // inside the cube
		Vector3d Xp = size/2. - X.abs();
		double d = std::min(Xp.x, std::min(Xp.y, Xp.z));

		return -1. * d;
	}

	double a = std::max(0., fabs(X.x) - size.x/2.);
	double b = std::max(0., fabs(X.y) - size.y/2.);
	double c = std::max(0., fabs(X.z) - size.z/2.);

	return sqrt(a*a + b*b +c*c);
}

Vector3d ParaxialBox::normal(const Vector3d& point) const
{
  Vector3d d = (point-corner).abs();
  Vector3d d2 = d + size;
  Vector3d n;
  double dmin = std::numeric_limits<double>::infinity();
  if (d.x < dmin)
  {
    dmin = d.x;
    n = Vector3d(-1,0,0);
  }
  if (d.y < dmin)
  {
    dmin = d.y;
    n = Vector3d(0,-1,0);
  }
  if (d.z < dmin)
  {
    dmin = d.z;
    n = Vector3d(0,0,-1);
  }
  if (d2.x < dmin)
  {
    dmin = d2.x;
    n = Vector3d(1,0,0);
  }
  if (d2.y < dmin)
  {
    dmin = d2.y;
    n = Vector3d(0,1,0);
  }
  if (d2.z < dmin)
  {
    dmin = d2.z;
    n = Vector3d(0,0,1);
  }

  return n;
}

std::string ParaxialBox::getDescription() const
{
	std::stringstream ss;
	ss << "ParaxialBox: " << std::endl
		 << "   corner: " << corner << std::endl
		 << "     size: " << size << std::endl;
	return ss.str();
};


} // namespace

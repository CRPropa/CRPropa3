#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

#include "crpropa/Candidate.h"
#include "crpropa/Vector3.h"

namespace crpropa
{
/**
 * \addtogroup Core
 * @{
 */

/**
 @class Surface
 @brief A geometrical surface

 Defines a surface. Can be queried if the candidate has crossed the surface in the last step.
 */

class Surface
{
	public:

	/**
		Returns the distance of a point to the surface. Negative on the one side,
		positive on the other. For closed surfaces it is negative on the inside.
	 */
    virtual double distance(const Vector3d& point) const = 0;
};

/**
 @class Plane
 @brief A plane given by a point x0 and two axes v1 and v2 with normal n = v1.cross(v2) or the normal n. Note that distance is negative on one side of the plane and positive on the other, depending on the orientation of the normal vector.
 */
class Plane: public Surface
{
	private:
		Vector3d x0, n;
	public:
		Plane(const Vector3d& _x0, const Vector3d& v1,const Vector3d& v2);
		Plane(const Vector3d& _x0, const Vector3d& _n);
    virtual double distance(const Vector3d &x) const;
};


/**
 @class Sphere
 @brief A sphere around point _center with radius _radius.
 */
class Sphere: public Surface
{
	private:
		Vector3d center;
		double radius;
	public:
		Sphere(const Vector3d& _center, double _radius);
    virtual double distance(const Vector3d &point) const;
};


///**
// @class Box
// @brief A box with arbitrary orientation and not necessarily perpendicular surfaces (Rhombohedron). For performance reasons use the ParaxialBox if applicable.
// */
//class Box: public Surface
//{
//	private:
//		Vector3d corner, x1, x2, x3, u, v, w;
//		std::vector<Plane> facets;
//		double A,B,C;
//	public:
//		Box(const Vector3d& _corner, const Vector3d& _x1,const Vector3d& _x2, const Vector3d& _x3);
//    virtual double distance(const Vector3d &point) const;
//};


/**
 @class ParaxialBox
 @brief A box with arbitrary orientation and not necessarily perpendicular surfaces (Rhombohedron)
 */
class ParaxialBox: public Surface
{
	private:
		Vector3d corner, size;
	public:
		ParaxialBox(const Vector3d& _corner, const Vector3d& _size);
    virtual double distance(const Vector3d &point) const;
};




/** @}*/
} // namespace crpropa

#endif // GEOMETRY_H

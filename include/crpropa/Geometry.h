#ifndef CRPROPA_GEOMETRY_H
#define CRPROPA_GEOMETRY_H

#include <string>

#include "crpropa/Candidate.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class Surface
 @brief A geometrical surface

 Defines a surface. Can be queried if the candidate has crossed the surface in the last step.
 */
class Surface  {
public:
	virtual ~Surface() = default;
	/** Returns the distance of a point to the surface. Negative on the one side,
 	 positive on the other. For closed surfaces it is negative on the inside.
	 @param point	vector corresponding to the point to which compute the distance
	 */
	virtual double distance(const Vector3d& point) const = 0;
	/** Returns the normal to the surface at a point. Negative on the one side,
	 positive on the other. For closed surfaces it is negative on the inside.
	 @param point	vector corresponding to the point to which compute the normal vector
	 */
	virtual Vector3d normal(const Vector3d& point) const = 0;
	virtual std::string getDescription() const {return "Surface without description.";};
};


/**
 @class Plane
 @brief A plane given by a point x0 and two axes v1 and v2 with normal n = v1.cross(v2) or the normal n. Note that distance is negative on one side of the plane and positive on the other, depending on the orientation of the normal vector.
 */
class Plane: public Surface {
private:
	Vector3d x0, n;
public:
	/** Constructor
	 * The plane is constructed from the given point x0 along the two axes v1 and v2
	 * @param x0  Point where to start construction
	 * @param v1  First axis of plane
	 * @param v2  Second axis of plane
	 */
	Plane(const Vector3d& x0, const Vector3d& v1,const Vector3d& v2);
	/** Constructor
	 * The plane is constucted from the given point x0 perpendicular on the normal n
	 * @param x0  Point where to start construction
	 * @param n  Normal
	 */
	Plane(const Vector3d& x0, const Vector3d& n);
	/**
	 * The distance is calculated with the normal n: return n.dot(x-x0);
	 */
	virtual double distance(const Vector3d &x) const;
	/// returns the normal of the plain independent of the point
	virtual Vector3d normal(const Vector3d& point=Vector3d(0,0,0)) const;
	virtual std::string getDescription() const;
};


/**
 @class Sphere
 @brief A sphere around point _center with radius _radius.
 */
class Sphere: public Surface {
private:
	Vector3d center;
	double radius;
public:
	/** Constructor
	 * @param center  Center point of the sphere
	 * @param radius  Radius of the sphere
	 */
	Sphere(const Vector3d& center, double radius);
	/**
	 * Distance to the sphere, negative on the inside
	 * @param point  Point to calculate the distance from the sphere from
	 */
	virtual double distance(const Vector3d &point) const;
	/** Returns the normal to the surface at a point. Negative on the inside,
	 positive on the other.
	 @param point  Vector corresponding to the point to which compute the normal vector
	 */
	virtual Vector3d normal(const Vector3d& point) const;
	virtual std::string getDescription() const;
};


/**
 @class ParaxialBox
 @brief A box with perpendicular surfaces aligned to the x,y,z-axes.
 */
class ParaxialBox: public Surface {
private:
	Vector3d corner, size;
public:
	/** Constructor
	 * @param corner  The Box will stretch beginning from corner in positive x,y,z direction
	 * @param size  The size of the box (all sides are equal)
	 */
	ParaxialBox(const Vector3d& corner, const Vector3d& size);
	/// Return distance to the Box for given point
	virtual double distance(const Vector3d &point) const;
	/// Return normal to the closest surface for given point
	virtual Vector3d normal(const Vector3d& point) const;
	virtual std::string getDescription() const;
};


/** @}*/
} // namespace crpropa

#endif // CRPROPA_GEOMETRY_H

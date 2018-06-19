#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>

#include "crpropa/Candidate.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"

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

class Surface : public Referenced
{
	public:

	/**
		Returns the distance of a point to the surface. Negative on the one side,
		positive on the other. For closed surfaces it is negative on the inside.
	 */
    virtual double distance(const Vector3d& point) const = 0;
	/**
		Returns the normal to the surface at a point. Negative on the one side,
		positive on the other. For closed surfaces it is negative on the inside.
	 */
    virtual Vector3d normal(const Vector3d& point) const = 0;
		virtual std::string getDescription() const {return "Surface without description.";};
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
    virtual Vector3d normal(const Vector3d& point) const;
		virtual std::string getDescription() const;
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
    virtual Vector3d normal(const Vector3d& point) const;
		virtual std::string getDescription() const;
};


/**
 @class ParaxialBox
 @brief A box with perpendicular surfaces aligned to the x,y,z-axes.
 */
class ParaxialBox: public Surface
{
	private:
		Vector3d corner, size;
	public:
		ParaxialBox(const Vector3d& _corner, const Vector3d& _size);
    virtual double distance(const Vector3d &point) const;
    virtual Vector3d normal(const Vector3d& point) const;
		virtual std::string getDescription() const;
};


/** @}*/
} // namespace crpropa

#endif // GEOMETRY_H

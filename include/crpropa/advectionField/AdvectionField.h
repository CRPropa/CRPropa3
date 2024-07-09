#ifndef CRPROPA_ADVECTIONFIELD_H
#define CRPROPA_ADVECTIONFIELD_H


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
	virtual Vector3d getField(const Vector3d &position, const double &time=0) const = 0;
	virtual double getDivergence(const Vector3d &position, const double &time=0) const = 0;
};


/**
 @class AdvectionFieldList
 @brief Advection field decorator implementing a superposition of fields.
 */
class AdvectionFieldList: public AdvectionField {
	std::vector<ref_ptr<AdvectionField> > fields;
public:
	void addField(ref_ptr<AdvectionField> field);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;
};


/**
 @class UniformAdvectionField
 @brief Advection field with one velocity/advection-field vector.
 */
class UniformAdvectionField: public AdvectionField {
	Vector3d value;
public:
	UniformAdvectionField(const Vector3d &value);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	std::string getDescription() const;
};


/**
@class ConstantSphericalAdvectionField
@brief Spherical advection field with a constant wind speed
*/

class ConstantSphericalAdvectionField: public AdvectionField {
	Vector3d origin; //origin of the advection sphere
	double vWind; // wind velocity
public:
	/** Constructor
	 @param origin	Origin of the advection field
	 @param vWind	Constant wind velocity

*/

	ConstantSphericalAdvectionField(const Vector3d origin, double vWind);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	void setOrigin(const Vector3d origin);
	void setVWind(double vMax);

	Vector3d getOrigin() const;
	double getVWind() const;

	std::string getDescription() const;

	
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
	/** Constructor
	@param origin 	Origin of the advection sphere
	@param radius 	Radius of the advection sphere
	@param vMax	Maximum wind velocity
	@param tau	Transition distance
	@param alpha	Tuning parameter
*/
	SphericalAdvectionField(const Vector3d origin, double radius, double vMax, double tau, double alpha);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	double getV(const double &r) const;

	void setOrigin(const Vector3d origin);
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

/**
 @class OneDimensionalCartesianShock
 @brief Advection field in x-direction with shock at x = 0 and width lShock approximated by tanh() 
		with variable compression ratio vUp/vDown
 */
class OneDimensionalCartesianShock: public AdvectionField {
	double compressionRatio; //compression ratio of shock
	double vUp; //upstream velocity 
	double lShock; //shock width
public:
/** Constructor
	@param compressionRatio //compression ratio of shock
	@param vUp //upstream velocity 
	@param lShock //shock width
*/
	OneDimensionalCartesianShock(double compressionRatio, double vUp, double lShock);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	void setComp(double compressionRatio);
	void setVup(double vUp);
	void setShockwidth(double lShock);

	double getComp() const;
	double getVup() const; 
	double getShockwidth() const;

	std::string getDescription() const;
};

/**
 @class OneDimensionalSphericalShock
 @brief Advection field in x-direction with shock at rShock and width lShock approximated by tanh() 
		with variable compression ratio ratio vUp/vDown
 */
class OneDimensionalSphericalShock: public AdvectionField {
	double compressionRatio; 	//compression ratio of shock
	double vUp; 	//upstream velocity 
	double lShock; 	//shock width
	double rShock; 	//shock radius
	bool coolUpstream; 	//flag for upstream cooling
public:
/** Constructor
	@param compressionRatio	//compression ratio of shock
	@param vUp 	//upstream velocity 
	@param lShock 	//shock width
	@param rShock 	//shock radius
	@param coolUpstream //flag for upstream cooling
*/
	OneDimensionalSphericalShock(double rShock, double vUp, double compressionRatio, double lShock, bool coolUpstream);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	void setComp(double compressionRatio);
	void setVup(double vUp);
	void setShockwidth(double lShock);
	void setShockRadius(double rShock);
	void setCooling(bool coolUpstream);

	double getComp() const; 
	double getVup() const;
	double getShockwidth() const;
	double getShockRadius() const;
	bool getCooling() const;

	std::string getDescription() const;
};

/**
 @class ObliqueAdvectionShock
 @brief Advection field in x-y-direction with shock at x = 0 and width x_sh approximated by tanh() 
		with variable compression ratio r_comp = vx_up/vx_down. The y component vy is not shocked 
		and remains constant. 
 */
class ObliqueAdvectionShock: public AdvectionField {
	double compressionRatio; //compression ratio of shock
	double vXUp; //upstream velocity x-component
	double vY; //constant velocity y-component
	double lShock; //shock width
	
public:
/** Constructor
	@param compressionRatio //compression ratio of shock
	@param vXUp //upstream velocity x-component
	@param vY //constant velocity y-component
	@param lShock //shock width
	
*/
	ObliqueAdvectionShock(double compressionRatio, double vXUp, double vY, double lShock);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	void setComp(double compressionRatio);
	void setVup(double vXUp);
	void setVy(double vY);
	void setShockwidth(double lShock);

	double getComp() const; 
	double getVup() const;
	double getVy() const;
	double getShockwidth() const;
	
	std::string getDescription() const;
};

/**
 @class SphericalAdvectionShock
 @brief Spherical advection with a constant velocity for r<r_0
	at the the shock the velocity drops to v_0/4. followed by
	a decrease proportional to 1/r^2.
*/

class SphericalAdvectionShock: public AdvectionField {
	Vector3d origin; // origin of the advection sphere
	double r_0; // position of the shock
	double v_0; // constant velocity
	double lambda; //transition width
	double r_rot; // normalization radius for rotation speed
	double v_phi; // rotation speed at r_rot

public:
	/** Constructor
	@param origin 	Origin of the advection sphere
	@param r_0 	Position of the shock
	@param v_0 	Constant velocity (r<<r_o)
	@param lambda 	Transition width / width of the shock
*/
	SphericalAdvectionShock(const Vector3d origin, double r_0, double v_0, double lambda);

	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	double g(double R) const;
	double g_prime(double R) const;

	void setOrigin(const Vector3d Origin);
	void setR0(double r);
	void setV0(double v);
	void setLambda(double l);
	void setRRot(double r);
	void setAzimuthalSpeed(double vPhi);

	Vector3d getOrigin() const;
	double getR0() const;
	double getV0() const;
	double getLambda() const;
	/**
	 * @param r Normalization radius for rotation speed
	*/	
	double getRRot() const;
	/**
	 * @param vPhi 	Rotation speed at r_rot
	*/	
	double getAzimuthalSpeed() const;

	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_ADVECTIONFIELD_H

#ifndef CRPROPA_TIMEDEPENDENTADVECTIONFIELD_H
#define CRPROPA_TIMEDEPENDENTADVECTIONFIELD_H

#include "crpropa/advectionField/AdvectionField.h"

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
 @class OneDimensionalTimeDependentShock
 @brief Advection field in x-direction with shock at x = x0 that propagates with a constant speed vsh = v_up - v_down
 		and width x_sh approximated by a tanh() with variable compression ratio r_comp = v_up/v_down.
        Pre- and postshock speeds as well as shock speed must be specified. 
        The shock position at t=0 can be specified, as well as the time t0 the shock starts to propagate.
 */
class OneDimensionalTimeDependentShock: public AdvectionField {
	double v_sh; // shock speed assuming xsh = vsh * t + xsh_0
	double v1; // speed behind the shock in lab frame
	double v0; // undisturbed speed in lab frame
	double l_sh; // shock width
    double x_sh0; // shock position at t = 0
	double t_sh0; // time the shock starts to propagate, before: vsh=0
public:/** Constructor
	@param v_sh; // shock speed
    @param v1; // speed behind the shock in lab frame
    @param v0; // undisturbed speed in lab frame
	@param l_sh; // shock width
*/
	OneDimensionalTimeDependentShock(double v_sh, double v1, double v0, double l_sh);

	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	void setShockSpeed(double v_sh);
	void setSpeeds(double v1, double v0);
	void setShockWidth(double l_sh);
	void setShockPosition(double x_sh0);
	void setShockTime(double t_sh0);

    double getVshock() const;
    double getV1() const;
    double getV0() const;
    double getShockWidth() const;
    double getShockPosition(double time) const;
    double getShockTime() const;
};

/**
 @class SedovTaylorBlastWave
 @brief Spherical advection field with shock at R(t), velocity vsh(t) and width l_sh approximated by tanh() 
		Wind solution is given by Kahn 1975, see also Drury 1983 for acceleration at the expanding shock
 */
class SedovTaylorBlastWave: public AdvectionField {
	double E0; 	    // energy of the explosion
	double rho0; 	// initial density
	double l_sh; 	// shock width
public:/** Constructor
	@param E0 // energy of the explosion
	@param rho0	// initial density
	
*/
	SedovTaylorBlastWave(double E0, double rho0, double l_sh);
	Vector3d getField(const Vector3d &position, const double &time=0) const;
	double getDivergence(const Vector3d &position, const double &time=0) const;

	void setShockWidth(double l_sh);
	void setEnergy(double E0);
	void setDensity(double rho0);

	double getShockRadius(double time) const;
    double getShockSpeed(double time) const;
    double getShockWidth() const;
    double getEnergy() const;
    double getDensity() const;
};

} // namespace crpropa

#endif // CRPROPA_TIMEDEPENDENTADVECTIONFIELD_H

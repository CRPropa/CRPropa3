#ifndef CRPROPA_JF12FIELDSOLENOIDAL_H
#define CRPROPA_JF12FIELDSOLENOIDAL_H

#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/Grid.h"
#include "crpropa/Units.h"

namespace crpropa {

/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class JF12FieldSolenoidal
 @brief JF12FieldSolenoidal galactic magnetic field model

 Implements a modified JF2012 Galactic magnetic field model.
 This implementation inherits most methods from the initial JF12Field,
 just the regular disk field and the poloidal halo "X field" are altered.

 A solenoidal transition for the disk field (which redirects the magnetic flux of the magnetic spiral field at the boundaries and lets the field strength tend to 0 in a continuous way)
 and an improved version of the JF12 X field (with parabolic field lines at z = 0 instead of sharp kinks)
 were added in order to smooth the field and avoid violations of magnetic flux conservation, see:
 Kleimann et al, 2018, arXiv:1809.07528, Solenoidal improvements for the JF12 Galactic magnetic field model.

 The turbulent field component is the exact same as in the initial JF12 implementation
 and should be used with care since the new regular and old turbulent field
 do not match in the transition regions of the spiral field.
 */

class JF12FieldSolenoidal: public JF12Field {
private:
	double zS; // height parameter of the modified X-field, field lines are parabolic for fabs(z) < zS
	double phi0Arms[11]; // azimuth angles in [-pi,pi] at which the dividing 8 spirals of the disk field intersect the (r1 = 5kpc)-ring at indices 1 to 8, remaining angles periodically filled.

	double phi0; // lower bound for azimuth phi integration to restore solenoidality of spiral field transition, arbitrary parameter

	// in order to restore solenoidality in the transition regions of the spiral field,
	// a phi-integral over the piecwise constant field strengths at r=5kpc has to be evaluated. Here,
	// we set H(phi) = phiCoeff[j] + bDisk[j] * phi as the result of this integration
	double phiCoeff[10]; // only 10 since the array holds the results for the integral from phi0Arms[0] to phi0Arms[1]; phi0Arms[0] to phi0Arms[2] etc.
	double corr; // correction term for enforcing H(phi0) = 0 afterwards which sets the lower boundary of the integration to phi0

	// inner and outer boundaries of disk field at r = 5 and 20 kpc
	double r1;
	double r2;

	// transitions region boundaries of disk field at r1 + delta, r2 - delta
	double r1s;
	double r2s;

public:
/** Constructor
	@param delta 	Transition width for the disk field such that the magnetic flux of the spiral field lines is redirected between r = 5 kpc and r = 5 kpc + delta as well as r = 20 kpc - delta and r = 20 kpc. The input parameter delta should be non-negative and smaller than 7.5 kpc. The default value is 3 kpc.
	@param zs 	Scale height for the X-shaped poloidal halo field such that the straight field lines are replaced by parabolas for abs(z) < zs. The input parameter zs should be non-negative. The default value is 500 pc.
*/
	JF12FieldSolenoidal(double delta = 3 * kpc, double zs = 0.5 * kpc);

	void setUseStriatedField(bool use); // override these for correct warning messages
	void setUseTurbulentField(bool use);

/** @brief Adjust the transition width of the disk field
	@param delta	The new transition width for the disk field with field strength transitions and flux redirection between r = 5 kpc and r = 5 kpc + delta as well as r = 20 kpc - delta and r = 20 kpc. Should be non-negative and smaller than 7.5 kpc.
	@return Void
*/
	void setDiskTransitionWidth(double delta);
	double getDiskTransitionWidth() const;


/** @brief Adjust the scale heigth of the X field
	@param zs	The new scale height of the X field with parabolic field lines for abs(z) < zs. Should be non-negative.
	@return Void
*/
	void setXScaleHeight(double zs);
	double getXScaleHeight() const;

	Vector3d getXField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const; // override old X and spiral field
	Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;

/** @brief Disable the transition of the spiral field strength to 0 at the outer boundary such that only the magnetic flux at the 5 kpc ring is redirected.
	Thus, the spiral field lines are continued to r = 20 kpc as in the initial JF12 field. You can reactivate the outer transition afterwards via setDiskTransitionWidth which sets both transition widths at the inner and outer boundary.
	@return Void
*/
	void deactivateOuterTransition();

/** @brief Evaluate the polynomial which provides the transition of the spiral field strength to zero in the transition regions.
	This transition is differentiable at the boundaries of the unaltered spiral field and contuinuous at the outer boundaries of the spiral field where the field strength goes to zero.
	@param r Distance of the current position to the z axis in the usual galactocentric cylindrical coordinates. Should be non-negative.
	@return The value of the transition polynomial at r if r is inside one of the transition regions. Otherwise, return (5 kpc)/r, i.e. the scaling of the spiral field, if r is inside the region where the spiral field remains unaltered.
*/
	double getDiskTransitionPolynomial(const double& r) const;

/** @brief Evaluate the derivative of the polynomial which provides the transition of the spiral field strength to zero in the transition regions. The derivative is needed to restore solenoidality in the transition region
	@param r 	Distance of the current position to the z axis in the usual galactocentric cylindrical coordinates. Should be non-negative.
	@return The value of the transition polynomial derivative at r if r is inside one of the transition regions. Otherwise, return 0.
*/
	double getDiskTransitionPolynomialDerivative(const double& r) const;

/** @brief Evaluate an angular azimuth integral over the piecewise constant spiral field strengths at r_1 = 5 kpc. This integral is needed to restore the solenoidality of the spiral field in the transition regions.
	@param r 	Distance of the current position to the z axis in the usual galactocentric cylindrical coordinates. Should be non-negative.
	@param phi	Azimuth angle of the current position in galactocentric cylindrical coordinates. Can be any double.
	@return The value of the azimuth integral over the spiral field strengths at r_1 = 5 kpc from at fixed phi0 to the current phi which is mapped back to r_1 = 5 kpc along the spiral field line passing through (r,phi).
*/
	double getHPhiIntegral(const double& r, const double& phi) const;

/** @brief Find the correct magnetic spiral arm for the current position and return its field strength at r_1 = 5 kpc
	@param r 	Distance of the current position to the z axis in the usual galactocentric cylindrical coordinates. Should be non-negative.
	@param phi	Azimuth angle of the current position in galactocentric cylindrical coordinates. Can be any double.
	@return The value of the spiral field strength b_j
*/
	double getSpiralFieldStrengthConstant(const double& r, const double& phi) const;
};
/** @} */

} // namespace crpropa

#endif // CRPROPA_JF12FIELDSOLENOIDAL_H

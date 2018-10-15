#ifndef CRPROPA_JF12FIELDSOLENOIDAL_H
#define CRPROPA_JF12FIELDSOLENOIDAL_H

#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/Grid.h"
#include <crpropa/Units.h>

namespace crpropa {

/**
 @class JF12FieldSolenoidal
 @brief JF12FieldSolenoidal galactic magnetic field model

 Implements a modified JF2012 magnetic field model, see documentation of JF12Field for basics.
 This implemementation inherits most methods and variables from the initial JF12Field, 
 just the regular disk and X field are altered.

 A solenoidal ring-spiral transition (of width delta = 3 kpc per default) for the disk field
 and a smooth X-field at z = 0 (parabolic field lines for abs(z) < zs; zs = 500 parsec per default)
 were added in order to smooth the field and avoid violations of magnetic flux conservation, see:
 Kleimann et al, 2018, arXiv:1809.07528, Solenoidal improvements for the JF12 Galactic magnetic field model.

 The turbulent field component is the exact same as in the initial JF12 implementation
 and should be used with care since the new regular and old turbulent field
 do not match in the ring-spiral transition region.
 */
 
class JF12FieldSolenoidal: public JF12Field {
private:

	// height parameter of the modified X-field, field lines are parabolic for fabs(z) < zS
	double zS;

	// angles at which the dividing 8 spirals of the disk field intersect the (r1 = 5kpc)-ring, 
	// longer array due to cyclic closure of the values
	double phi0Arms[11];

	// lower bound for phi integration, somewhat arbitrary parameter
	double phi0;
	// corresponding index in phi0Arms such that phi0Arms[idx0-1] < phi0 < phi0Arms[idx0]
	int idx0;

	// phi integration H(phi) = phiCoeff[j] + bDisk[j] * phi
	double phiCoeff[10]; // for phi0Arms[j-1] < phi < phi0Arms[j], hence only 10
	double corr; // correction term enforcing H(phi0) = 0

	// inner and outer boundaries of disk field at r = 5 and 20 kpc
	double r1;
	double r2;
	
	// transitions region boundaries of disk field at r1 + delta, r2 - delta
	double r1s;
	double r2s;

	// initial JF12 field strengths of the 8 spiral arms  (index 1 to 8) at r = 5 kpc,
	// cyclically closed s.t. b[0] = b[8], b[9] = b[1], b[10] = b[2]
	double bDiskCyclicClosure[11]; 

public:
	JF12FieldSolenoidal(double delta= 3 * kpc, double zs= 0.5 * kpc);
	
	void setUseStriated(bool use); // override these for correct warning messages
	void setUseTurbulent(bool use);

	// set and get transition width at inner and outer boundary of spiral field
	void setDelta(double d);
	double getDelta() const;

	// set and get scale heigth for parabolic X field lines
	void setZs(double z);
	double getZs() const;

	// override old X and spiral field
	Vector3d getXField(const double& r, const double& z, const double& sinPhi, const double& cosPhi) const;
	Vector3d getDiskField(const double& r, const double& z, const double& phi, const double& sinPhi, const double& cosPhi) const;

	// continue spiral field to r=20kpc without transition at the outer boundary
	// one may reactivate the transition afterwards with setDelta()
	void deactivateOuterTransition();

	// transition polynomial p_delta(r) for the spiral field
	double p(const double& r) const;

	// transition polynomial derivative
	double q(const double& r) const;

	// evaluate H-Integral needed for solenoidality of the spiral field
	double PhiIntegralH(const double& r, const double& phi) const;

	// return field strength b_j of initial field at r = 5 kpc for current spiral arm
	double getSpiralStrength(const double& r, const double& phi) const;
};

} // namespace crpropa

#endif // CRPROPA_JF12FIELDSOLENOIDAL_H

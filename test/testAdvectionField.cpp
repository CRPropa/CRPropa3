#include "crpropa/advectionField/AdvectionField.h"
#include "crpropa/advectionField/TimeDependentAdvectionField.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include "gtest/gtest.h"
#include <stdexcept>
#include <cmath>

namespace crpropa {

TEST(testUniformAdvectionField, SimpleTest) {
	UniformAdvectionField A(Vector3d(-1, 5, 3));
	Vector3d a = A.getField(Vector3d(1, 0, 0));
	double D = A.getDivergence(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(a.x, -1);
	EXPECT_DOUBLE_EQ(a.y, 5);
	EXPECT_DOUBLE_EQ(a.z, 3);
	EXPECT_DOUBLE_EQ(D, 0.);
}

TEST(testAdvectionFieldList, SimpleTest) {
	// Test a list of three advection fields
	AdvectionFieldList A;
	A.addField(new UniformAdvectionField(Vector3d(1, 0, 0)));
	A.addField(new UniformAdvectionField(Vector3d(0, 2, 0)));
	A.addField(new UniformAdvectionField(Vector3d(0, 0, 3)));
	Vector3d a = A.getField(Vector3d(0.));
	double D = A.getDivergence(Vector3d(1, 2, 3));
	EXPECT_DOUBLE_EQ(a.x, 1);
	EXPECT_DOUBLE_EQ(a.y, 2);
	EXPECT_DOUBLE_EQ(a.z, 3);
	EXPECT_DOUBLE_EQ(D, 0.);
}

TEST(testConstantSphericalAdvectionField, SimpleTest) {
	
	Vector3d origin(1, 0, 0);
	double V_wind(10);
	
	ConstantSphericalAdvectionField A(origin, V_wind);
	
	// Check the properties of the advection field
	EXPECT_DOUBLE_EQ(A.getOrigin().x, origin.x);
	EXPECT_DOUBLE_EQ(A.getOrigin().y, origin.y);
	EXPECT_DOUBLE_EQ(A.getOrigin().z, origin.z);

	EXPECT_DOUBLE_EQ(A.getVWind(), V_wind);

	// Field should be radial with center (1,0,0)
	Vector3d Pos(2, 1, 1);
	Vector3d V = A.getField(Pos);

	EXPECT_DOUBLE_EQ(V.x, 10.*pow(3, -0.5));
	EXPECT_DOUBLE_EQ(V.y, 10.*pow(3, -0.5));
	EXPECT_DOUBLE_EQ(V.z, 10.*pow(3, -0.5));
	
	// Divergence should be 2*V/r
	Vector3d Pos2(2, 0, 0);
	double D = A.getDivergence(Pos2);
	EXPECT_DOUBLE_EQ(D, 2*10);
}

TEST(testSphericalAdvectionField, SimpleTest) {

	Vector3d origin(1, 0, 0);
	double R_max(10);
	double V_max(1000);
	double tau(3.);
	double alpha(2.);

	SphericalAdvectionField A(origin, R_max, V_max, tau, alpha);

	// Check the properties of the advection field
	EXPECT_DOUBLE_EQ(A.getOrigin().x, origin.x);
	EXPECT_DOUBLE_EQ(A.getOrigin().y, origin.y);
	EXPECT_DOUBLE_EQ(A.getOrigin().z, origin.z);
	
	EXPECT_DOUBLE_EQ(A.getRadius(), R_max);
	EXPECT_DOUBLE_EQ(A.getVMax(), V_max);
	EXPECT_DOUBLE_EQ(A.getTau(), tau);
	EXPECT_DOUBLE_EQ(A.getAlpha(), alpha);

	// Field/divergence should be zero for R>10
	EXPECT_DOUBLE_EQ(A.getField(Vector3d(12,0,0)).getR(), 0.);
	EXPECT_DOUBLE_EQ(A.getDivergence(Vector3d(12,0,0)), 0.);

	// Check field and divergence
	Vector3d Pos(2, 0, 0);
	Vector3d a = A.getField(Pos);
	Vector3d a0 = a.getUnitVector();
	double d = A.getDivergence(Pos);

	EXPECT_DOUBLE_EQ(a0.x, 1.);
	EXPECT_DOUBLE_EQ(a0.y, 0.);
	EXPECT_DOUBLE_EQ(a0.z, 0.);
	EXPECT_DOUBLE_EQ(a.getR(), V_max*(1.-exp(-1./tau)) );

	EXPECT_NEAR(d, 1044.624919, 1e-5);

	// Check asymptotic of the Field
	EXPECT_NEAR(A.getField(Vector3d(11, 0, 0)).getR(), 1000., 1e-4);
	
	// Divergence should be 2*V_max/r for r>>1
	EXPECT_NEAR(A.getDivergence(Vector3d(11, 0, 0)), 2*1000./10., 1e-4);
}


TEST(testSphericalAdvectionShock, SimpleTest) {

	Vector3d origin(0, 0, 0);
	double R_0(10);
	double V_0(1000);
	double lambda(0.1);
	double R_rot(1.);
	double V_rot(100);
	

	SphericalAdvectionShock A(origin, R_0, V_0, lambda);

	// Check the properties of the advection field
	EXPECT_DOUBLE_EQ(A.getOrigin().x, origin.x);
	EXPECT_DOUBLE_EQ(A.getOrigin().y, origin.y);
	EXPECT_DOUBLE_EQ(A.getOrigin().z, origin.z);
	
	EXPECT_DOUBLE_EQ(A.getR0(), R_0);
	EXPECT_DOUBLE_EQ(A.getV0(), V_0);
	EXPECT_DOUBLE_EQ(A.getLambda(), lambda);

	// Default values for R_Rot is R_0
	// Default value for Azimuthal speed is 0
	EXPECT_DOUBLE_EQ(A.getRRot(), R_0);
	EXPECT_DOUBLE_EQ(A.getAzimuthalSpeed(), 0.);	

	// Field should drop to 0.625*V_0 at the shock 
	// That's a difference to the analytic shock model where it should
	// drop to v_r(R_0)=0.25*V_0.
	EXPECT_DOUBLE_EQ(A.getField(Vector3d(R_0,0,0)).getR(), 0.625*V_0);

	// Field divergence should be zero for R>>R_0
	EXPECT_NEAR(A.getDivergence(Vector3d(15,0,0)), 0., 1e-10);

	// Check field and divergence
	Vector3d Pos(2, 0, 0);
	Vector3d a = A.getField(Pos);
	Vector3d a0 = a.getUnitVector();
	double d = A.getDivergence(Pos);

	EXPECT_DOUBLE_EQ(a0.x, 1.);
	EXPECT_DOUBLE_EQ(a0.y, 0.);
	EXPECT_DOUBLE_EQ(a0.z, 0.);
	EXPECT_DOUBLE_EQ(a.getR(), V_0);
	
	//Set explicitely the azimuthal rotation speed 
	A.setRRot(R_rot);
	A.setAzimuthalSpeed(V_rot);
	
	EXPECT_DOUBLE_EQ(A.getRRot(), R_rot);
	EXPECT_DOUBLE_EQ(A.getAzimuthalSpeed(), V_rot);	
	
	Vector3d pos(1., 0., 0.);
	Vector3d f = A.getField(pos);
	EXPECT_DOUBLE_EQ(f.x, 1000.);
	EXPECT_DOUBLE_EQ(f.y, 100.);
	EXPECT_DOUBLE_EQ(f.z, 0.);	

	
}

TEST(testOneDimensionalTimeDependentShock, SimpleTest) {

    double V_sh(1);
    double V_1(3./4.);
    double V_0(0.);
    double L_sh(.01);
    double X_sh0(0);
    double T_sh0(0);

    OneDimensionalTimeDependentShock A(V_sh, V_1, V_0, L_sh);

    // Check the properties of the advection field
    EXPECT_DOUBLE_EQ(A.getVshock(), V_sh);
    EXPECT_DOUBLE_EQ(A.getV1(), V_1);
    EXPECT_DOUBLE_EQ(A.getV0(), V_0);
    EXPECT_DOUBLE_EQ(A.getShockWidth(), L_sh);
    EXPECT_DOUBLE_EQ(A.getShockPosition(0.), X_sh0);
    EXPECT_DOUBLE_EQ(A.getShockTime(), T_sh0);

    // Field should be zero for t=0
    EXPECT_DOUBLE_EQ(A.getField(Vector3d(1,0,0)).getR(), 0.);
    EXPECT_DOUBLE_EQ(A.getDivergence(Vector3d(1,0,0)), 0.);

    // Shock position at t=1
    double xsh = A.getShockPosition(10.);
    EXPECT_DOUBLE_EQ(xsh, V_sh * 10.);

    // Preshock and postshock speeds:
    EXPECT_DOUBLE_EQ(A.getField(Vector3d(xsh+5.,0,0), 10.).getR(), V_0);
    EXPECT_DOUBLE_EQ(A.getField(Vector3d(xsh-5.,0,0), 10.).getR(), V_1);

}

TEST(testSedovTaylorBlastWave, SimpleTest) {

    double E0(1);
    double rho0(1);
    double L_sh(0.01);

    SedovTaylorBlastWave A(E0, rho0, L_sh);

    // Check the properties of the advection field
    EXPECT_DOUBLE_EQ(A.getEnergy(), E0);
    EXPECT_DOUBLE_EQ(A.getDensity(), rho0);
    EXPECT_DOUBLE_EQ(A.getShockWidth(), L_sh);

    // Field should be zero for t=0
    EXPECT_DOUBLE_EQ(A.getField(Vector3d(1,0,0)).getR(), 0.);
    EXPECT_DOUBLE_EQ(A.getDivergence(Vector3d(1,0,0)), 0.);

    // Shock position at t=1
    double R = A.getShockRadius(1.);
    EXPECT_DOUBLE_EQ(R, pow(E0 / rho0, 1./5.));

    // Shock speed at t=1
    double V = A.getShockSpeed(1.);
    EXPECT_DOUBLE_EQ(V, 2. / 5. * pow(E0 / rho0, 1. / 5.));

    // Check field 
    Vector3d Pos(2, 0, 0);
    Vector3d a = A.getField(Pos);

    EXPECT_DOUBLE_EQ(a.getR(), 0.5 * 2. * pow(2., 8) * 2. / 5. * pow(E0 / rho0, 1. / 5.) * (1 - tanh( (2 - 1) * pow(E0 / rho0, 1. / 5.) / L_sh) ));

    // Check asymptotic of the Field
    EXPECT_NEAR(A.getField(Vector3d(11, 0, 0)).getR(), 0, 1e-4);
    EXPECT_NEAR(A.getDivergence(Vector3d(11, 0, 0)), 0, 1e-4);

}

} //namespace crpropa

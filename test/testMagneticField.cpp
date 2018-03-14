#include <stdexcept>

#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testUniformMagneticField, SimpleTest) {
	UniformMagneticField B(Vector3d(-1, 5, 3));
	Vector3d b = B.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(b.getX(), -1);
	EXPECT_DOUBLE_EQ(b.getY(), 5);
	EXPECT_DOUBLE_EQ(b.getZ(), 3);
}

TEST(testMagneticDipoleField, SimpleTest) {
	// Test magnetic dipole
	// mu0 / (4*M_PI) * m / r^3 (2*cos(theta)*e_r + sin(theta)*e_theta)
	MagneticDipoleField B(Vector3d(0,0,0), Vector3d(0,0,1), 1);
	Vector3d b1 = B.getField(Vector3d(0, 0, 1)); // theta = 0
	Vector3d b2 = B.getField(Vector3d(1, 0, 0)); // theta = 0
	EXPECT_NEAR(b1.getX(), 0, 1E-8);
	EXPECT_NEAR(b1.getY(), 0, 1E-8);
	EXPECT_NEAR(b1.getZ(), mu0 / (4*M_PI) * 2, 1E-8);
	EXPECT_NEAR(b2.getX(), 0, 1E-8);
	EXPECT_NEAR(b2.getY(), 0, 1E-8);
	EXPECT_NEAR(b2.getZ(), -1 * mu0 / (4*M_PI), 1E-8);
}

#ifdef CRPROPA_HAVE_MUPARSER
TEST(testRenormalizeMagneticField, simpleTest) {
	ref_ptr<UniformMagneticField> field = new UniformMagneticField(Vector3d(2*nG, 0, 0));
	RenormalizeMagneticField modField(field, "B^2-1*nG");
	Vector3d b = modField.getField(Vector3d(5, 5, 5));
	EXPECT_NEAR(b.getR(), 3*nG, 0.001);
}
#endif

TEST(testMagneticFieldList, SimpleTest) {
	// Test a list of three magnetic fields
	MagneticFieldList B;
	B.addField(new UniformMagneticField(Vector3d(1, 0, 0)));
	B.addField(new UniformMagneticField(Vector3d(0, 2, 0)));
	B.addField(new UniformMagneticField(Vector3d(0, 0, 3)));
	Vector3d b = B.getField(Vector3d(0.));
	EXPECT_DOUBLE_EQ(b.x, 1);
	EXPECT_DOUBLE_EQ(b.y, 2);
	EXPECT_DOUBLE_EQ(b.z, 3);
}

TEST(testMagneticFieldEvolution, SimpleTest) {
	// Test if this decorator scales the underlying field as (1+z)^m
	ref_ptr<UniformMagneticField> B = new UniformMagneticField(Vector3d(1,0,0));
	double z = 1.2;
	double m = 3;
	MagneticFieldEvolution Bz(B, m);

	// scaled field
	Vector3d b = Bz.getField(Vector3d(0,0,0), z);
	EXPECT_DOUBLE_EQ(b.x, pow(1+z, m));

	// unscaled field
	b = Bz.getField(Vector3d(0,0,0));
	EXPECT_DOUBLE_EQ(b.x, 1);
}

#ifdef CRPROPA_HAVE_FFTW3F
TEST(testVectorFieldGrid, Turbulence_bmean_brms) {
	// Test for zero mean: <B> = 0
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;

	ref_ptr<VectorGrid> grid = new VectorGrid(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(grid, Brms, lMin, lMax);

	double precision = 1e-7;
	Vector3f bMean = meanFieldVector(grid);
	EXPECT_NEAR(0, bMean.x, precision);
	EXPECT_NEAR(0, bMean.y, precision);
	EXPECT_NEAR(0, bMean.z, precision);
	EXPECT_NEAR(1, rmsFieldStrength(grid), precision);
}

TEST(testVectorFieldGrid, Turbulence_seed) {
	// Test if seeding produces 2 identical fields
	size_t n = 64;
	double spacing = 1 * Mpc;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	double index = -11. / 3.;
	int seed = 753;

	ref_ptr<VectorGrid> grid1 = new VectorGrid(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(grid1, Brms, lMin, lMax, index, seed);

	ref_ptr<VectorGrid> grid2 = new VectorGrid(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(grid2, Brms, lMin, lMax, index, seed);

	Vector3d pos(22 * Mpc);
	EXPECT_FLOAT_EQ(grid1->interpolate(pos).x, grid2->interpolate(pos).x);
}

TEST(testVectorFieldGrid, turbulence_Exceptions) {
	// Test exceptions
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double brms = 1;
	ref_ptr<VectorGrid> grid = new VectorGrid(Vector3d(0, 0, 0), n, spacing);

	// should be fine
	EXPECT_NO_THROW(initTurbulence(grid, brms, 2 * spacing, 8 * spacing));
	// lMin too small
	EXPECT_THROW(initTurbulence(grid, brms, 1.5 * spacing, 8 * spacing),
			std::runtime_error);
	// lMin > lMax
	EXPECT_THROW(initTurbulence(grid, brms, 8.1 * spacing, 8 * spacing),
			std::runtime_error);
	// lMax too large
	EXPECT_THROW(initTurbulence(grid, brms, 2 * spacing, 33 * spacing),
			std::runtime_error);
}
#endif // CRPROPA_HAVE_FFTW3F

class EchoMagneticField: public MagneticField {
public:
	Vector3d getField(const Vector3d &position) const {
		return position;
	}
};

TEST(testPeriodicMagneticField, Exceptions) {
	ref_ptr<EchoMagneticField> f = new EchoMagneticField();
	ref_ptr<PeriodicMagneticField> p = new PeriodicMagneticField(f,
			Vector3d(10000, 10000, 10000), Vector3d(1000, 1000, 1000), true);

	// box 0, 0, 0
	Vector3d v = p->getField(Vector3d(1000, 2000, 3000));
	EXPECT_DOUBLE_EQ(0, v.x);
	EXPECT_DOUBLE_EQ(1000, v.y);
	EXPECT_DOUBLE_EQ(2000, v.z);

	// box 1, 2, 3
	v = p->getField(Vector3d(12000, 23000, 34000));
	EXPECT_DOUBLE_EQ(9000, v.x);
	EXPECT_DOUBLE_EQ(2000, v.y);
	EXPECT_DOUBLE_EQ(7000, v.z);

	// box -1, -2, -3
	v = p->getField(Vector3d(0, -10000, -20000));
	EXPECT_DOUBLE_EQ(1000, v.x);
	EXPECT_DOUBLE_EQ(9000, v.y);
	EXPECT_DOUBLE_EQ(1000, v.z);

}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

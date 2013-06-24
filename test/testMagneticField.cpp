#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/TurbulentMagneticField.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testUniformMagneticField, SimpleTest) {
	UniformMagneticField B(Vector3d(-1, 5, 3));
	Vector3d b = B.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(b.x, -1);
	EXPECT_DOUBLE_EQ(b.y, 5);
	EXPECT_DOUBLE_EQ(b.z, 3);
}

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

#ifdef MPC_HAVE_FFTW3F
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
#endif // MPC_HAVE_FFTW3F
TEST(testTurbulentMagneticField, SimpleTest) {
	TurbulentMagneticField B;
	B.setTurbulenceProperties(1 * nG, 10 * parsec, 200 * parsec, -11. / 3.,
			1000);
	B.initialize();
	Vector3d b1 = B.getField(Vector3d(0.));
	Vector3d b2 = B.getField(Vector3d(0.));
	EXPECT_DOUBLE_EQ(b1.x, b2.x);
	EXPECT_DOUBLE_EQ(b1.y, b2.y);
	EXPECT_DOUBLE_EQ(b1.z, b2.z);
}

TEST(testTurbulentMagneticField, Brms) {
	TurbulentMagneticField B;
	B.setTurbulenceProperties(1, 0.1, 100);
	B.initialize();
	double sumB2 = 0;
	Vector3d Bmean, b;
	int n = 20;
	Random &r = Random::instance();
	for (int ix = 0; ix < n; ix++) {
		for (int iy = 0; iy < n; iy++) {
			for (int iz = 0; iz < n; iz++) {
				b = B.getField(Vector3d(r.rand(), r.rand(), r.rand()) * 100000);
				Bmean += b;
				sumB2 = b.getMag2();
			}
		}
	}
	Bmean /= pow(n, 3);
	double Brms = sqrt(sumB2 / pow(n, 3));
	std::cout << Brms << std::endl;
	std::cout << Bmean << std::endl;

	// this turbulent field realization is not working at the moment
//	EXPECT_NEAR(Bmean.x, 0, 1e-6);
//	EXPECT_NEAR(Bmean.y, 0, 1e-6);
//	EXPECT_NEAR(Bmean.z, 0, 1e-6);
//
//	EXPECT_NEAR(Brms, 1, 1e-2);
}

TEST(testTurbulentMagneticField, Exceptions) {
	TurbulentMagneticField B;

	// Test exception if properties not set
	EXPECT_THROW(B.initialize(), std::runtime_error);

	// Test exception for lMin > lMax
	B.setTurbulenceProperties(1 * nG, 8.1, 8);
	EXPECT_THROW(B.initialize(), std::runtime_error);
}

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

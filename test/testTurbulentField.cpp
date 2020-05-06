#include <stdexcept>

#include "crpropa/Grid.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/GridTools.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/turbulentField/GridTurbulence.h"
#include "crpropa/magneticField/turbulentField/PlaneWaveTurbulence.h"
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testTurbulentField, correlationLength) {
	double l_bo = 100;
	auto tf = TurbulentField(1*muG, 5./3., 4., l_bo);
	auto Lc = tf.getCorrelationLength();
    EXPECT_NEAR(Lc, 0.498*l_bo, 0.001*l_bo);
}

TEST(testPlaneWaveTurbulence, correlationLength) {
    double Brms = 1*muG;
    double lMin = 1*kpc;
    double lMax = 800*kpc;
	double l_bo = 100*kpc;
    double s = 5/3.;
	double q = 4.;
    
    auto tf = PlaneWaveTurbulence(Brms, s, q, l_bo, lMin, lMax);
    auto Lc = tf.getCorrelationLength();
    EXPECT_NEAR(Lc, 0.498*l_bo, 1*kpc);
}

#ifdef CRPROPA_HAVE_FFTW3F

TEST(testSimpleGridTurbulence, oldFunctionForCrrelationLength) { //TODO: remove in future
	double lMin = 1*kpc;
	double lMax = 1*Gpc;
	double alpha = -11/3.;
	auto Lc = turbulentCorrelationLength(lMin, lMax, alpha);
    EXPECT_NEAR(Lc, lMax/5, 1*Mpc);
}

TEST(testVectorFieldGrid, Turbulence_bmean_brms) {
	// Test for zero mean: <B> = 0
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	double sindex = 5/3.;

	ref_ptr<Grid3f> grid = new Grid3f(Vector3d(0, 0, 0), n, spacing);
    auto tf = SimpleGridTurbulence(grid, Brms, sindex, lMin, lMax);

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
	double sindex = 5/3.;
	int seed = 753;

	ref_ptr<Grid3f> grid1 = new Grid3f(Vector3d(0, 0, 0), n, spacing);
    auto tf1 = SimpleGridTurbulence(grid1, Brms, sindex, lMin, lMax, seed);

	ref_ptr<Grid3f> grid2 = new Grid3f(Vector3d(0, 0, 0), n, spacing);
    auto tf2 = SimpleGridTurbulence(grid2, Brms, sindex, lMin, lMax, seed);

	Vector3d pos(22 * Mpc);
	EXPECT_FLOAT_EQ(grid1->interpolate(pos).x, grid2->interpolate(pos).x);
}

TEST(testVectorFieldGrid, turbulence_Exceptions) {
	// Test exceptions
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double brms = 1;
	ref_ptr<Grid3f> grid = new Grid3f(Vector3d(0, 0, 0), n, spacing);

	// should be fine
	EXPECT_NO_THROW(initTurbulence(grid, brms, 2 * spacing, 8 * spacing));
	// lMin too small
	EXPECT_THROW(initTurbulence(grid, brms, 1.5 * spacing, 8 * spacing),
			std::runtime_error);
	// lMin > lMax
	EXPECT_THROW(initTurbulence(grid, brms, 8.1 * spacing, 8 * spacing),
			std::runtime_error);
	// lMax too large
	EXPECT_THROW(initTurbulence(grid, brms, 2 * spacing, 65 * spacing),
			std::runtime_error);
}

TEST(testGridTurbulence, Turbulence_seed) {
	// Test if seeding produces 2 identical fields
	size_t n = 64;
	double spacing = 1 * Mpc;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	double s = 5/3.;
	double q = 4;
	double l_bo = lMax/6;
	int seed = 137;

	ref_ptr<Grid3f> grid1 = new Grid3f(Vector3d(0, 0, 0), n, spacing);
    auto tf1 = GridTurbulence(grid1, Brms, s, q, l_bo, lMin, lMax, seed);

	ref_ptr<Grid3f> grid2 = new Grid3f(Vector3d(0, 0, 0), n, spacing);
    auto tf2 = GridTurbulence(grid2, Brms, s, q, l_bo, lMin, lMax, seed);

	Vector3d pos(22 * Mpc);
	EXPECT_FLOAT_EQ(tf1.getField(pos).x, tf2.getField(pos).x);
}
#endif // CRPROPA_HAVE_FFTW3F

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

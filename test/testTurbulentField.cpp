#include <stdexcept>

#include "crpropa/Grid.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/GridTools.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testTurbulentField, correlationLength) {
	double l_bo = 100;
	auto tf = TurbulentField(1*muG, 5./3., 4., l_bo);
	auto Lc = tf.getCorrelationLength();
    EXPECT_NEAR(Lc, 0.498*l_bo, 0.001*l_bo);
}

#ifdef CRPROPA_HAVE_FFTW3F

TEST(testSimpleGridTurbulence, correlationLength) {
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

	ref_ptr<Grid3f> grid = new Grid3f(Vector3d(0, 0, 0), n, spacing);
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

	ref_ptr<Grid3f> grid1 = new Grid3f(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(grid1, Brms, lMin, lMax, index, seed);

	ref_ptr<Grid3f> grid2 = new Grid3f(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(grid2, Brms, lMin, lMax, index, seed);

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
#endif // CRPROPA_HAVE_FFTW3F

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

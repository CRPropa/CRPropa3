#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/magneticField/PeriodicGrid.h"
#include "mpc/magneticField/TurbulentMagneticField.h"
#include "mpc/Units.h"
#include "mpc/Common.h"

#include "gtest/gtest.h"
#include <stdexcept>

namespace mpc {

TEST(testUniformMagneticField, SimpleTest) {
	UniformMagneticField B(Vector3d(-1, 5, 3));
	Vector3d b = B.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(b.x, -1);
	EXPECT_DOUBLE_EQ(b.y, 5);
	EXPECT_DOUBLE_EQ(b.z, 3);
}

TEST(testVectorFieldGrid, PeriodicClamp) {
	// Test correct determination of lower and upper neighbor
	int lo, hi;

	periodicClamp(23.12, 8, lo, hi);
	EXPECT_EQ(7, lo);
	EXPECT_EQ(0, hi);

	periodicClamp(-23.12, 8, lo, hi);
	EXPECT_EQ(0, lo);
	EXPECT_EQ(1, hi);
}

TEST(testVectorFieldGrid, SimpleTest) {
	// Test construction and parameters
	VectorFieldGrid B(Vector3d(1., 2., 3.), 4, 2.5);
	EXPECT_TRUE(Vector3d(1., 2., 3.) == B.getOrigin());
	EXPECT_EQ(4, B.getNx());
	EXPECT_EQ(4, B.getNy());
	EXPECT_EQ(4, B.getNz());
	EXPECT_DOUBLE_EQ(2.5, B.getSpacing());
}

TEST(testVectorFieldGrid, Interpolation) {
	// Explicitly test trilinear interpolation
	double spacing = 2.793;
	VectorFieldGrid grid(Vector3d(0.), 3, spacing);
	grid.get(0, 0, 1) = Vector3f(1.7, 0., 0.); // set one value

	Vector3d b;
	b = grid.interpolate(Vector3d(0, 0, 1) * spacing);
	EXPECT_FLOAT_EQ(1.7, b.x);

	b = grid.interpolate(Vector3d(0, 0, 0.9) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9, b.x);

	b = grid.interpolate(Vector3d(0, 0, 1.1) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9, b.x);

	b = grid.interpolate(Vector3d(0, 0.15, 0.9) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.85, b.x);

	b = grid.interpolate(Vector3d(0, 2.15, 0.9) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.15, b.x);
}

TEST(testVectorFieldGrid, Scale) {
	// Test scaling a field
	VectorFieldGrid grid(Vector3d(0.), 3, 1);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid.get(ix, iy, iz) = Vector3f(1, 0, 0);

	scale(&grid, 5);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				EXPECT_FLOAT_EQ(5, grid.interpolate(Vector3d(0.7, 0, 0.1)).x);
}

TEST(testVectorFieldGrid, Periodicity) {
	// Test for periodic boundaries: B(x+a*n) = B(x)
	size_t n = 3;
	double spacing = 3;
	double size = n * spacing;
	VectorFieldGrid grid(Vector3d(0.), n, spacing);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid.get(ix, iy, iz) = Vector3f(iz + ix, iy * iz, ix - iz * iy);

	Vector3d pos(1.2, 2.3, 0.7);
	Vector3f b = grid.interpolate(pos);
	Vector3f b2 = grid.interpolate(pos + Vector3d(1, 0, 0) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	b2 = grid.interpolate(pos + Vector3d(0, 5, 0) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	b2 = grid.interpolate(pos + Vector3d(0, 0, -2) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);
}

TEST(testVectorFieldGrid, DumpLoad) {
	// Dump and load a field grid
	VectorFieldGrid B1(Vector3d(0.), 3, 1);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				B1.get(ix, iy, iz) = Vector3f(1, 2, 3);
	dump(&B1, getDataPath("../test/testDump.raw"));

	VectorFieldGrid B2(Vector3d(0.), 3, 1);
	load(&B2, getDataPath("../test/testDump.raw"));

	for (int ix = 0; ix < 3; ix++) {
		for (int iy = 0; iy < 3; iy++) {
			for (int iz = 0; iz < 3; iz++) {
				Vector3f b1 = B1.get(ix, iy, iz);
				Vector3f b2 = B2.get(ix, iy, iz);
				EXPECT_FLOAT_EQ(b1.x, b2.x);
				EXPECT_FLOAT_EQ(b1.y, b2.y);
				EXPECT_FLOAT_EQ(b1.z, b2.z);
			}
		}
	}
}

TEST(testVectorFieldGrid, DumpLoadTxt) {
	// Dump and load a field grid
	VectorFieldGrid B1(Vector3d(0.), 3, 1);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				B1.get(ix, iy, iz) = Vector3f(ix, iy, iz) * nG;
	dumpTxt(&B1, getDataPath("../test/testDump.txt"), 1e4);

	VectorFieldGrid B2(Vector3d(0.), 3, 1);
	loadTxt(&B2, getDataPath("../test/testDump.txt"), 1e-4);

	for (int ix = 0; ix < 3; ix++) {
		for (int iy = 0; iy < 3; iy++) {
			for (int iz = 0; iz < 3; iz++) {
				Vector3f b1 = B1.get(ix, iy, iz);
				Vector3f b2 = B2.get(ix, iy, iz);
				EXPECT_FLOAT_EQ(b1.x, b2.x);
				EXPECT_FLOAT_EQ(b1.y, b2.y);
				EXPECT_FLOAT_EQ(b1.z, b2.z);
			}
		}
	}
}

TEST(testVectorFieldGrid, Speed) {
	// Dump and load a field grid
	VectorFieldGrid grid(Vector3d(0.), 3, 3);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid.get(ix, iy, iz) = Vector3f(1, 2, 3);

	Vector3d b;
	for (int i = 0; i < 100000; i++)
		b = grid.interpolate(Vector3d(i));
}

#ifdef MPC_HAVE_FFTW3F
TEST(testVectorFieldGrid, Turbulence_bmean_brms) {
	// Test for zero mean: <B> = 0
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	VectorFieldGrid grid(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(&grid, Brms, lMin, lMax);

	double precision = 1e-7;
	Vector3f bMean = meanFieldStrength(&grid);
	EXPECT_NEAR(0, bMean.x, precision);
	EXPECT_NEAR(0, bMean.y, precision);
	EXPECT_NEAR(0, bMean.z, precision);
	EXPECT_NEAR(1, rmsFieldStrength(&grid), precision);
}

TEST(testMagneticFieldGrid, Turbulence_seed) {
	// Test if seeding produces 2 identical fields
	size_t n = 64;
	double spacing = 1 * Mpc;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	double index = -11. / 3.;
	int seed = 753;

	VectorFieldGrid grid1(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(&grid1, Brms, lMin, lMax, index, seed);

	VectorFieldGrid grid2(Vector3d(0, 0, 0), n, spacing);
	initTurbulence(&grid2, Brms, lMin, lMax, index, seed);

	Vector3d pos(22 * Mpc);
	EXPECT_FLOAT_EQ(grid1.interpolate(pos).x, grid2.interpolate(pos).x);
}

TEST(testMagneticFieldGrid, turbulence_Exceptions) {
	// Test exceptions
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double brms = 1;
	VectorFieldGrid grid(Vector3d(0, 0, 0), n, spacing);

	// should be fine
	EXPECT_NO_THROW(initTurbulence(&grid, brms, 2 * spacing, 8 * spacing));
	// lMin too small
	EXPECT_THROW(initTurbulence(&grid, brms, 1.5 * spacing, 8 * spacing),
			std::runtime_error);
	// lMin > lMax
	EXPECT_THROW(initTurbulence(&grid, brms, 8.1 * spacing, 8 * spacing),
			std::runtime_error);
	// lMax too large
	EXPECT_THROW(initTurbulence(&grid, brms, 2 * spacing, 33 * spacing),
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
	Random random;
	for (int ix = 0; ix < n; ix++) {
		for (int iy = 0; iy < n; iy++) {
			for (int iz = 0; iz < n; iz++) {
				b = B.getField(
						Vector3d(random.rand(), random.rand(), random.rand())
								* 100000);
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

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

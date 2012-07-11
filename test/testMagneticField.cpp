#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/magneticField/TurbulentMagneticFieldGrid.h"
#include "mpc/magneticField/TurbulentMagneticField.h"
#include "mpc/Units.h"

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

TEST(testMagneticFieldGrid, PeriodicClamp) {
	// Test correct determination of lower and upper neighbor
	int lo, hi;

	periodicClamp(23.12, 8, lo, hi);
	EXPECT_EQ(7, lo);
	EXPECT_EQ(0, hi);

	periodicClamp(-23.12, 8, lo, hi);
	EXPECT_EQ(0, lo);
	EXPECT_EQ(1, hi);
}

TEST(testMagneticFieldGrid, SimpleTest) {
	// Test construction and parameters
	MagneticFieldGrid B(Vector3d(1., 2., 3.), 10, 4);
	EXPECT_TRUE(Vector3d(1., 2., 3.) == B.getGridOrigin());
	EXPECT_EQ(4, B.getGridSamples());
	EXPECT_DOUBLE_EQ(10, B.getGridSize());
	EXPECT_DOUBLE_EQ(2.5, B.getGridSpacing());
}

TEST(testMagneticFieldGrid, Interpolation) {
	// Explicitly test trilinear interpolation
	MagneticFieldGrid B(Vector3d(0.), 9, 3);
	B.get(0, 0, 1) = Vector3f(1.7, 0., 0.);

	double spacing = B.getGridSpacing();
	EXPECT_FLOAT_EQ(1.7, B.getField(Vector3d(0, 0, 1) * spacing).x);

	EXPECT_FLOAT_EQ(1.7 * 0.9, B.getField(Vector3d(0, 0, 0.9) * spacing).x);
	EXPECT_FLOAT_EQ(1.7 * 0.9, B.getField(Vector3d(0, 0, 1.1) * spacing).x);

	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.85,
			B.getField(Vector3d(0, 0.15, 0.9) * spacing).x);
	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.15,
			B.getField(Vector3d(0, 2.15, 0.9) * spacing).x);
}

TEST(testMagneticFieldGrid, Normalization) {
	// If the field is uniform, (interpolated) field values should be the same everywhere
	MagneticFieldGrid B(Vector3d(0.), 3, 3);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				B.get(ix, iy, iz) = Vector3f(1, 0, 0);

	EXPECT_FLOAT_EQ(1, B.getField(Vector3d(0.7, 0, 0.1)).x);
	EXPECT_FLOAT_EQ(1, B.getField(Vector3d(0, 1.3, 0.9)).x);

	B.normalize(2);

	EXPECT_FLOAT_EQ(2, B.getField(Vector3d(0.7, 0, 0.1)).x);
	EXPECT_FLOAT_EQ(2, B.getField(Vector3d(0, 1.3, 0.9)).x);
}

TEST(testMagneticFieldGrid, Periodicity) {
	// Test for periodic boundaries: B(x+a*n) = B(x)
	MagneticFieldGrid B(Vector3d(0.), 9, 3);

	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				B.get(ix, iy, iz) = Vector3f(iz + ix, iy * iz, ix - iz * iy);

	double size = B.getGridSize();
	Vector3d pos(1.2, 2.3, 0.7);
	Vector3f b = B.getField(pos);
	Vector3f b2 = B.getField(pos + Vector3d(1, 0, 0) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	b2 = B.getField(pos + Vector3d(0, 5, 0) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	b2 = B.getField(pos + Vector3d(0, 0, -2) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);
}

TEST(testMagneticFieldGrid, DumpLoad) {
	// Dump and load a field grid
	MagneticFieldGrid B1(Vector3d(0.), 3, 3);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				B1.get(ix, iy, iz) = Vector3f(1, 2, 3);
	B1.dump("test/testDump.raw");

	MagneticFieldGrid B2(Vector3d(0.), 3, 3);
	B2.load("test/testDump.raw");

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

TEST(testMagneticFieldGrid, Speed) {
	// Dump and load a field grid
	MagneticFieldGrid B(Vector3d(0.), 3, 3);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				B.get(ix, iy, iz) = Vector3f(1, 2, 3);

	Vector3d b;
	Vector3d pos(1.,1.,1.);
	for (int i = 0; i < 1000000; i++)
		b = B.getField(pos);
}

TEST(testTurbulentMagneticFieldGrid, Bmean) {
	// Test for zero mean: <B> = 0
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	TurbulentMagneticFieldGrid B(Vector3d(0, 0, 0), 10 * Mpc, n, 2 * spacing,
			8 * spacing, -11. / 3., 1.);

	Vector3d b(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++)
				b += B.getField(Vector3d(ix, iy, iz) * spacing);
	b /= n * n * n;

	double precision = 1e-9;
	EXPECT_NEAR(b.x, 0, precision);
	EXPECT_NEAR(b.y, 0, precision);
	EXPECT_NEAR(b.z, 0, precision);
}

TEST(testTurbulentMagneticFieldGrid, Brms) {
	// Test for correct RMS strength: <B^2> = Brms^2
	size_t n = 100;
	TurbulentMagneticFieldGrid B(Vector3d(0, 0, 0), 10 * Mpc, n);
	double spacing = B.getGridSpacing();
	B.setTurbulenceProperties(2 * spacing, 8 * spacing, -11. / 3.);
	B.initialize();
	B.normalize(1. / B.getRMSFieldStrength());
	double brms = 0;
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++)
				brms += B.getField(Vector3d(ix, iy, iz) * spacing).getMag2();
	brms = sqrt(brms / n / n / n);

	double precision = 1e-7;
	EXPECT_NEAR(brms, 1, precision);
}

TEST(testTurbulentMagneticFieldGrid, Exceptions) {
	// Test exceptions
	TurbulentMagneticFieldGrid B(Vector3d(0, 0, 0), 10, 64);
	double spacing = B.getGridSpacing();

	// should be fine
	B.setTurbulenceProperties(2 * spacing, 8 * spacing, -11. / 3.);
	EXPECT_NO_THROW(B.initialize());

	// lMin too small
	B.setTurbulenceProperties(1.5 * spacing, 8 * spacing, -11. / 3.);
	EXPECT_THROW(B.initialize(), std::runtime_error);

	// lMin > lMax
	B.setTurbulenceProperties(8.1 * spacing, 8 * spacing, -11. / 3.);
	EXPECT_THROW(B.initialize(), std::runtime_error);

	// lMax too large
	B.setTurbulenceProperties(2 * spacing, 33 * spacing, -11. / 3.);
	EXPECT_THROW(B.initialize(), std::runtime_error);
}

TEST(testTurbulentMagneticField, SimpleTest) {
	TurbulentMagneticField B;
	B.setTurbulenceProperties(1 * nG, 10 * parsec, 200 * parsec, -11./3., 1000);
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
				b = B.getField(Vector3d(random.rand(), random.rand(), random.rand()) * 100000);
				Bmean += b;
				sumB2 = b.getMag2();
			}
		}
	}
	Bmean /= pow(n, 3);
	double Brms = sqrt(sumB2 / pow(n, 3));
	std::cout << Brms << std::endl;
	std::cout << Bmean << std::endl;

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

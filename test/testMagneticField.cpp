#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/magneticField/TurbulentMagneticField.h"
#include "mpc/Units.h"

#ifdef MPC_HAVE_GADGET
#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#endif

#include "gtest/gtest.h"

namespace mpc {

TEST(testUniformMagneticField, SimpleTest) {
	UniformMagneticField B(Vector3d(-1, 5, 3));
	Vector3d b = B.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(b.x, -1);
	EXPECT_DOUBLE_EQ(b.y, 5);
	EXPECT_DOUBLE_EQ(b.z, 3);
}

TEST(testTurbulentMagneticFieldGrid, PeriodicBoundaries) {
	// Test turbulent field for periodic boundaries: B(x+a*n) = B(x)
	size_t n = 64;
	TurbulentMagneticField bField(Vector3d(0, 0, 0), n, 1);
	bField.initialize(2, 8, 1, -11. / 3.);

	Vector3d pos(1.1, 2.1, 3.1);
	Vector3d b = bField.getField(pos);
	Vector3d b1 = bField.getField(pos + Vector3d(n, 0, 0));
	Vector3d b2 = bField.getField(pos + Vector3d(0, n, 0));
	Vector3d b3 = bField.getField(pos + Vector3d(0, 0, n));

	EXPECT_FLOAT_EQ(b.x, b1.x);
	EXPECT_FLOAT_EQ(b.y, b1.y);
	EXPECT_FLOAT_EQ(b.z, b1.z);

	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	EXPECT_FLOAT_EQ(b.x, b3.x);
	EXPECT_FLOAT_EQ(b.y, b3.y);
	EXPECT_FLOAT_EQ(b.z, b3.z);
}

TEST(testTurbulentMagneticFieldGrid, ZeroMean) {
	// Test turbulent field for zero mean: <B> = 0
	size_t n = 64;
	TurbulentMagneticField bField(Vector3d(0, 0, 0), n, 1);
	bField.initialize(2, 8, 1, -11. / 3.);

	Vector3d b(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++)
				b += bField.getField(Vector3d(ix, iy, iz));

	b /= n * n * n;
	double precision = 1e-7;
	EXPECT_NEAR(b.x, 0, precision);
	EXPECT_NEAR(b.y, 0, precision);
	EXPECT_NEAR(b.z, 0, precision);
}

TEST(testTurbulentMagneticFieldGrid, Brms) {
	// Test turbulent field for correct RMS strength: <B^2> = Brms^2
	size_t n = 64;
	TurbulentMagneticField bField(Vector3d(0, 0, 0), n, 1);
	bField.initialize(2, 8, 1, -11. / 3.);

	double brms = 0;
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++)
				brms += bField.getField(Vector3d(ix, iy, iz)).getMag2();
	brms = sqrt(brms / n / n / n);

	double precision = 1e-7;
	EXPECT_NEAR(brms, 1, precision);
}

TEST(testSPHTurbulentMagneticField, construction) {
	SPHTurbulentMagneticField bField(Vector3d(80, 80, 80) * Mpc, 64, 40./64);
	bField.initialize(2, 8, 1, -11. / 3.);
	bField.modulate("/home/walz/software/mpc/data/mhd1.db");
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

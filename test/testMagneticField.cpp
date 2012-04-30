#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/magneticField/magneticFieldGrid.h"
#include "mpc/magneticField/turbulentMagneticFieldGrid.h"
#include "mpc/Units.h"

#include "gtest/gtest.h"
#include "fftw3.h"
#include <iostream>

namespace mpc {

TEST(testUniformMagneticField, SimpleTest) {
	UniformMagneticField B(Vector3(-1, 5, 3));
	Vector3 b = B.getField(Vector3(1, 0, 0));
	EXPECT_DOUBLE_EQ(b.x(), -1);
	EXPECT_DOUBLE_EQ(b.y(), 5);
	EXPECT_DOUBLE_EQ(b.z(), 3);
}

TEST(testTurbulentMagneticFieldGrid, PeriodicBoundaries) {
	// B(x+a*n) = B(x)

	size_t n = 64;
	TurbulentMagneticFieldGrid B(Vector3(0, 0, 0), n, 1, 2, 8, 1, -11. / 3.);

	Vector3 pos(1.1, 2.1, 3.1);
	Vector3 b = B.getField(pos);
	Vector3 b1 = B.getField(pos + Vector3(n, 0, 0));
	Vector3 b2 = B.getField(pos + Vector3(0, n, 0));
	Vector3 b3 = B.getField(pos + Vector3(0, 0, n));

	EXPECT_FLOAT_EQ(b.x(), b1.x());
	EXPECT_FLOAT_EQ(b.y(), b1.y());
	EXPECT_FLOAT_EQ(b.z(), b1.z());

	EXPECT_FLOAT_EQ(b.x(), b2.x());
	EXPECT_FLOAT_EQ(b.y(), b2.y());
	EXPECT_FLOAT_EQ(b.z(), b2.z());

	EXPECT_FLOAT_EQ(b.x(), b3.x());
	EXPECT_FLOAT_EQ(b.y(), b3.y());
	EXPECT_FLOAT_EQ(b.z(), b3.z());
}

TEST(testTurbulentMagneticFieldGrid, ZeroMean) {
	// <B> = 0

	size_t n = 64;
	TurbulentMagneticFieldGrid B(Vector3(0, 0, 0), n, 1, 2, 8, 1, -11. / 3.);

	Vector3 b(0, 0, 0);
	for (unsigned int ix = 0; ix < n; ix++)
		for (unsigned int iy = 0; iy < n; iy++)
			for (unsigned int iz = 0; iz < n; iz++)
				b += B.getField(Vector3(ix, iy, iz));
	b /= n * n * n;

	double precision = 1e-12;
	EXPECT_NEAR(b.x(), 0, precision);
	EXPECT_NEAR(b.y(), 0, precision);
	EXPECT_NEAR(b.z(), 0, precision);
}

TEST(testTurbulentMagneticFieldGrid, Brms) {
	// <B^2> = Brms^2

	size_t n = 64;
	TurbulentMagneticFieldGrid B(Vector3(0, 0, 0), n, 1*Mpc, 2*Mpc, 8*Mpc, 1., -11. / 3.);

	double brms = 0;
	for (unsigned int ix = 0; ix < n; ix++)
		for (unsigned int iy = 0; iy < n; iy++)
			for (unsigned int iz = 0; iz < n; iz++)
				brms += B.getField(Vector3(ix, iy, iz)*Mpc).mag2();
	brms = sqrt(brms / n / n / n);

	double precision = 1e-12;
	EXPECT_NEAR(brms, 1, precision);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

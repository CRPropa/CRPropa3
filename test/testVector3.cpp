#include "crpropa/Vector3.h"
#include "gtest/gtest.h"

namespace crpropa {

TEST(Vector3, division) {
	Vector3d v(10);
	v /= 10;
	EXPECT_DOUBLE_EQ(v.x, 1);
	EXPECT_DOUBLE_EQ(v.y, 1);
	EXPECT_DOUBLE_EQ(v.z, 1);
	v = Vector3d(10) / Vector3d(2); // element-wise division
	EXPECT_DOUBLE_EQ(v.x, 5);
	EXPECT_DOUBLE_EQ(v.y, 5);
	EXPECT_DOUBLE_EQ(v.z, 5);
}

TEST(Vector3, mod) {
	Vector3d v(10.1, 10.2, 10.3);
	v %= 10.2;
	EXPECT_NEAR(v.x, 10.1, 1e-10); // mod doesn't preserve double precision
	EXPECT_NEAR(v.y, 0, 1e-10);
	EXPECT_NEAR(v.z, 0.1, 1e-10);
}

TEST(Vector3, dot) {
	double a = Vector3d(1, 0, 0).dot(Vector3d(0, 1, 0));
	EXPECT_DOUBLE_EQ(a, 0);
	double b = Vector3d(-1, 10, 2).dot(Vector3d(5, 1, -3));
	EXPECT_DOUBLE_EQ(b, -1);
}

TEST(Vector3, cross) {
	Vector3d v = Vector3d(1, 0, 0).cross(Vector3d(0, 1, 0));
	EXPECT_DOUBLE_EQ(v.x, 0);
	EXPECT_DOUBLE_EQ(v.y, 0);
	EXPECT_DOUBLE_EQ(v.z, 1);
}

TEST(Vector3, angle) {
	double a = Vector3d(1, 1, 0).getAngleTo(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(a, 45 * M_PI / 180);
	double b = Vector3d(0, 0, 1).getAngleTo(Vector3d(0, 0, 1));
	EXPECT_DOUBLE_EQ(b, 0);
}

TEST(Vector3, magnitude) {
	Vector3d v = Vector3d(1, 2, -2);
	EXPECT_DOUBLE_EQ(v.getMag(), 3);
	EXPECT_DOUBLE_EQ(v.getMag2(), 9);
}

TEST(Vector3, distance) {
	double a = Vector3d(10, 0, 10).getDistanceTo(Vector3d(10, 0, 0));
	EXPECT_DOUBLE_EQ(a, 10);
}

TEST(Vector3, rotation) {
	Vector3d v(10, 0, 0);
	v.rotate(Vector3d(0, 2, 0), M_PI);
	EXPECT_NEAR(v.x, -10, 1e-9);
	// rotation doesn't preserve double precision
	EXPECT_NEAR(v.y, 0, 1e-9);
	EXPECT_NEAR(v.z, 0, 1e-9);

	v = Vector3d(10, 0, 0);
	v.rotate(Vector3d(1, 0, 0), M_PI); // no rotation
	EXPECT_NEAR(v.x, 10, 1e-9);
	EXPECT_NEAR(v.y, 0, 1e-9);
	EXPECT_NEAR(v.z, 0, 1e-9);

	v = Vector3d(5, 2, 7);
	v.rotate(Vector3d(1, 8, -4), 2 * M_PI); // full rotation
	EXPECT_NEAR(v.x, 5, 1e-9);
	EXPECT_NEAR(v.y, 2, 1e-9);
	EXPECT_NEAR(v.z, 7, 1e-9);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

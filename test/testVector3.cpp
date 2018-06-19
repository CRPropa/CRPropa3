#include "crpropa/Vector3.h"
#include "gtest/gtest.h"

namespace crpropa {

TEST(Vector3, comparison) {
	EXPECT_TRUE(Vector3d(1, 2, 3) == Vector3d(1, 2, 3));
	EXPECT_FALSE(Vector3d(1, 2, 3) == Vector3d(2, 3, 4));
}

TEST(Vector3, multiplication) {
	Vector3d v(1);
	v *= 10;
	EXPECT_DOUBLE_EQ(v.x, 10);
	EXPECT_DOUBLE_EQ(v.y, 10);
	EXPECT_DOUBLE_EQ(v.z, 10);
	v = Vector3d(1) * Vector3d(2); // element-wise multiplication
	EXPECT_DOUBLE_EQ(v.x, 2);
	EXPECT_DOUBLE_EQ(v.y, 2);
	EXPECT_DOUBLE_EQ(v.z, 2);
}

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
	EXPECT_TRUE(v == Vector3d(0, 0, 1));
}

TEST(Vector3, angle) {
	// 45 degrees
	double a = Vector3d(1, 1, 0).getAngleTo(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(a, 45 * M_PI / 180);
	// perpendicular vectors
	double b = Vector3d(0, 0, 1).getAngleTo(Vector3d(0, 0, 1));
	EXPECT_DOUBLE_EQ(b, 0);
}

TEST(Vector3, unitVectors) {
	Vector3d v = Vector3d(2, 0, 0);
	Vector3d er = v.getUnitVector();
	Vector3d et = v.getUnitVectorTheta();
	Vector3d ep = v.getUnitVectorPhi();

	// trigonometrical functions don't preserve double precision
	double eps = 1e-16;
	EXPECT_NEAR(er.x, 1, eps);
	EXPECT_NEAR(er.y, 0, eps);
	EXPECT_NEAR(er.z, 0, eps);

	EXPECT_NEAR(et.x, 0, eps);
	EXPECT_NEAR(et.y, 0, eps);
	EXPECT_NEAR(et.z, -1, eps);

	EXPECT_NEAR(ep.x, 0, eps);
	EXPECT_NEAR(ep.y, 1, eps);
	EXPECT_NEAR(ep.z, 0, eps);
}

TEST(Vector3, magnitude) {
	Vector3d v = Vector3d(1, 2, -2);
	EXPECT_DOUBLE_EQ(v.getR(), 3);
	EXPECT_DOUBLE_EQ(v.getR2(), 9);
}

TEST(Vector3, distance) {
	double a = Vector3d(10, 0, 10).getDistanceTo(Vector3d(10, 0, 0));
	EXPECT_DOUBLE_EQ(a, 10);
}

TEST(Vector3, rotation) {
	Vector3d v, w;

	// rotation by 180 degrees
	v = Vector3d(10, 0, 0);
	w = v.getRotated(Vector3d(0, 2, 0), M_PI);
	EXPECT_NEAR(w.x, -10, 1e-9);
	EXPECT_NEAR(w.y, 0, 1e-9);
	EXPECT_NEAR(w.z, 0, 1e-9);

	// rotation axis parallel to vector --> no rotation
	v = Vector3d(10, 0, 0);
	w = v.getRotated(Vector3d(1, 0, 0), M_PI);
	EXPECT_NEAR(w.x, 10, 1e-9);
	EXPECT_NEAR(w.y, 0, 1e-9);
	EXPECT_NEAR(w.z, 0, 1e-9);

	// rotation by 360 degrees
	v = Vector3d(5, 2, 7);
	w = v.getRotated(Vector3d(1, 8, -4), 2 * M_PI);
	EXPECT_NEAR(w.x, 5, 1e-9);
	EXPECT_NEAR(w.y, 2, 1e-9);
	EXPECT_NEAR(w.z, 7, 1e-9);

	// rotation by zero degrees
	v = Vector3d(1, 0, 0);
	w = v.getRotated(Vector3d(0,1,0), 0);

	EXPECT_NEAR(w.x, 1, 1e-9);
	EXPECT_NEAR(w.y, 0, 1e-9);
	EXPECT_NEAR(w.z, 0, 1e-9);

	// rotation around zero axis 
	v = Vector3d(1, 0, 0);
	Vector3d a = v.cross(v);
	w = v.getRotated(a, 0);

	EXPECT_NEAR(w.x, 1, 1e-9);
	EXPECT_NEAR(w.y, 0, 1e-9);
	EXPECT_NEAR(w.z, 0, 1e-9);




}

TEST(Vector3, parallelPerpendicular) {
	Vector3d v(3, 2, 1);
	Vector3d w(0, 1, 0);

	Vector3d vpara = v.getParallelTo(w);
	Vector3d vperp = v.getPerpendicularTo(w);

	EXPECT_DOUBLE_EQ(vpara.x, 0);
	EXPECT_DOUBLE_EQ(vpara.y, 2);
	EXPECT_DOUBLE_EQ(vpara.z, 0);
	EXPECT_DOUBLE_EQ(vperp.x, 3);
	EXPECT_DOUBLE_EQ(vperp.y, 0);
	EXPECT_DOUBLE_EQ(vperp.z, 1);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

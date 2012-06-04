#include "mpc/module/BreakCondition.h"
#include "mpc/Candidate.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(MinimumEnergy, Continue) {
	MinimumEnergy minEnergy(5);
	Candidate candidate;
	candidate.current.setEnergy(5.1);
	minEnergy.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(MinimumEnergy, Stop) {
	MinimumEnergy minEnergy(5);
	Candidate candidate;
	candidate.current.setEnergy(4.9);
	minEnergy.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(MaximumTrajectoryLength, Continue) {
	MaximumTrajectoryLength maxLength(10);
	Candidate candidate;
	candidate.setTrajectoryLength(9.9);
	maxLength.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(MaximumTrajectoryLength, stop) {
	MaximumTrajectoryLength maxLength(10);
	Candidate candidate;
	candidate.setTrajectoryLength(10.1);
	maxLength.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(SmallObserverSphere, Continue) {
	SmallObserverSphere obs(Vector3d(0, 0, 0), 1);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(2, 0, 0));
	obs.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(SmallObserverSphere, detect) {
	SmallObserverSphere obs(Vector3d(0, 0, 0), 1);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(0.1, 0.5, -0.1));
	obs.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
	EXPECT_TRUE(candidate.hasProperty("Detected"));
}

TEST(SmallObserverSphere, limitStep) {
	SmallObserverSphere obs(Vector3d(0, 0, 0), 1);
	Candidate candidate;
	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3d(0, 0, 2));
	obs.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1);
}

TEST(CubicBoundary, inside) {
	CubicBoundary cube(Vector3d(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(9, 5, 5));
	cube.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(CubicBoundary, outside) {
	CubicBoundary cube(Vector3d(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(10.1, 5, 5));

	cube.setMakeInactive(false);
	cube.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
	EXPECT_TRUE(candidate.hasProperty("OutOfBounds"));

	cube.setMakeInactive(true);
	cube.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(CubicBoundary, limitStep) {
	CubicBoundary cube(Vector3d(0, 0, 0), 10);
	cube.setLimitStep(true, 1);
	Candidate candidate;
	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3d(5, 5, 0.5));
	cube.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST(SphericalBoundary, inside) {
	SphericalBoundary sphere(Vector3d(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(9, 0, 0));
	sphere.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
	EXPECT_FALSE(candidate.hasProperty("OutOfBounds"));
}

TEST(SphericalBoundary, outside) {
	SphericalBoundary sphere(Vector3d(0, 0, 0), 10, "PassedGalacticBorder");

	Candidate candidate;
	candidate.current.setPosition(Vector3d(0, -10.1, 0));

	sphere.setMakeInactive(false);
	sphere.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
	EXPECT_TRUE(candidate.hasProperty("PassedGalacticBorder"));

	sphere.setMakeInactive(true);
	sphere.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(SphericalBoundary, limitStep) {
	SphericalBoundary sphere(Vector3d(0, 0, 0), 10);
	sphere.setLimitStep(true, 1);
	Candidate candidate;
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3d(0, 0, 9.5));
	sphere.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST(EllipsoidalBoundary, inside) {
	EllipsoidalBoundary ellipsoid(Vector3d(-5, 0, 0), Vector3d(5, 0, 0), 15);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(3, 2, 0));
	ellipsoid.process(&candidate);
	EXPECT_FALSE(candidate.hasProperty("OutOfBounds"));
}

TEST(EllipsoidalBoundary, outside) {
	EllipsoidalBoundary ellipsoid(Vector3d(-5, 0, 0), Vector3d(5, 0, 0), 15);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(0, 25, 0));

	ellipsoid.setMakeInactive(false);
	ellipsoid.process(&candidate);
	EXPECT_TRUE(candidate.hasProperty("OutOfBounds"));
	EXPECT_TRUE(candidate.isActive());

	ellipsoid.setMakeInactive(true);
	ellipsoid.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(EllipsoidalBoundary, limitStep) {
	EllipsoidalBoundary ellipsoid(Vector3d(-5, 0, 0), Vector3d(5, 0, 0), 15);
	ellipsoid.setLimitStep(true, 0.5);
	Candidate candidate;
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3d(7, 0, 0));
	ellipsoid.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

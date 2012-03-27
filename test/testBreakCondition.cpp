#include "mpc/module/BreakCondition.h"
#include "mpc/Candidate.h"
#include "gtest/gtest.h"

namespace mpc {

class ModuleTest: public testing::Test {
protected:
};

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

TEST(MaximumTrajectoryLength, Stop) {
	MaximumTrajectoryLength maxLength(10);
	Candidate candidate;
	candidate.setTrajectoryLength(10.1);
	maxLength.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(SmallObserverSphere, Continue) {
	SmallObserverSphere obs(Vector3(0, 0, 0), 1);
	Candidate candidate;
	candidate.current.setPosition(Vector3(2, 0, 0));
	obs.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(SmallObserverSphere, Detect) {
	SmallObserverSphere obs(Vector3(0, 0, 0), 1);
	Candidate candidate;
	candidate.current.setPosition(Vector3(0.1, 0.5, -0.1));
	obs.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
	EXPECT_TRUE(candidate.hasProperty("Detected"));
}

TEST(SmallObserverSphere, LimitStep) {
	SmallObserverSphere obs(Vector3(0, 0, 0), 1);
	Candidate candidate;
	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(0, 0, 2));
	obs.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1);
}

TEST(CubicBoundary, Inside) {
	CubicBoundary cube(Vector3(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3(9, 5, 5));
	cube.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(CubicBoundary, Outside) {
	CubicBoundary cube(Vector3(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3(10.1, 5, 5));
	cube.process(&candidate);
	EXPECT_FALSE(candidate.isActive());

	candidate.current.setPosition(Vector3(5, -0.1, 5));
	candidate.setActive(true);
	cube.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(CubicBoundary, LimitStep) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, 1);
	Candidate candidate;
	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(5, 5, 0.5));
	cube.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST(SphericalBoundary, Inside) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3(9, 0, 0));
	sphere.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(SphericalBoundary, Outside) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10);
	Candidate candidate;
	candidate.current.setPosition(Vector3(0, -10.1, 0));
	sphere.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(SphericalBoundary, LimitStep) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, 1);
	Candidate candidate;
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3(0, 0, 9.5));
	sphere.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST(EllipsoidalBoundary, Inside) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15);
	Candidate candidate;
	candidate.current.setPosition(Vector3(3, 2, 0));
	ellipsoid.process(&candidate);
	EXPECT_TRUE(candidate.isActive());
}

TEST(EllipsoidalBoundary, Outside) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15);
	Candidate candidate;
	candidate.current.setPosition(Vector3(0, 25, 0));
	ellipsoid.process(&candidate);
	EXPECT_FALSE(candidate.isActive());
}

TEST(EllipsoidalBoundary, LimitStep) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15, 0.5);
	Candidate candidate;
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3(7, 0, 0));
	ellipsoid.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

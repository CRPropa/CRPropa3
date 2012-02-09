#include "mpc/module/BreakCondition.h"
#include "mpc/Candidate.h"
#include "gtest/gtest.h"

namespace mpc {

TEST(testMinimumEnergy, Continue) {
	MinimumEnergy minEnergy(5);
	Candidate candidate;
	candidate.current.setEnergy(5.1);
	minEnergy.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST(testMinimumEnergy, Break) {
	MinimumEnergy minEnergy(5);
	Candidate candidate;
	candidate.current.setEnergy(4.9);
	minEnergy.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testMaximumTrajectoryLength, Continue) {
	MaximumTrajectoryLength maxLength(10);
	Candidate candidate;
	candidate.setTrajectoryLength(9.9);
	maxLength.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST(testMaximumTrajectoryLength, Break) {
	MaximumTrajectoryLength maxLength(10);
	Candidate candidate;
	candidate.setTrajectoryLength(10.1);
	maxLength.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testBreakCondition, SmallObserverSphere) {
	SmallObserverSphere obs(Vector3(0, 0, 0), 1);
	Candidate candidate;

	candidate.current.setPosition(Vector3(2, 0, 0));
	obs.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);

	candidate.setStatus(Candidate::Active);
	candidate.current.setPosition(Vector3(0.1, 0.5, -0.1));
	obs.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Detected);

	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(0, 0, 2));
	obs.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1);
}

TEST(testCubicBoundary, Inside) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(9, 5, 5));
	cube.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST(testCubicBoundary, Above) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(10.1, 5, 5));
	cube.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testCubicBoundary, Below) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(5, -0.1, 5));
	cube.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testCubicBoundary, LimitStep) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, 1, Candidate::Stopped);
	Candidate candidate;
	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(5, 5, 0.5));
	cube.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST(testSphericalBoundary, Inside) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(9, 0, 0));
	sphere.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST(testSphericalBoundary, Outside) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(0, -10.1, 0));
	sphere.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testSphericalBoundary, LimitStep) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, 1, Candidate::Stopped);
	Candidate candidate;
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3(0, 0, 9.5));
	sphere.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST(testEllipsoidalBoundary, Inside) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(3, 2, 0));
	ellipsoid.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST(testEllipsoidalBoundary, Outside) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15, Candidate::Stopped);
	Candidate candidate;
	candidate.current.setPosition(Vector3(0, 25, 0));
	ellipsoid.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testEllipsoidalBoundary, LimitStep) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15, 0.5, Candidate::Stopped);
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

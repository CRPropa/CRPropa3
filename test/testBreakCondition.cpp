#include "mpc/module/BreakCondition.h"
#include "mpc/Candidate.h"
#include "gtest/gtest.h"

namespace mpc {

class ModuleTest: public testing::Test {
protected:
	Candidate candidate;
	std::vector<Candidate *> secondaries;
};

TEST_F(ModuleTest, Continue) {
	MinimumEnergy minEnergy(5);
	candidate.current.setEnergy(5.1);
	minEnergy.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST_F(ModuleTest, Break) {
	MinimumEnergy minEnergy(5);
	candidate.current.setEnergy(4.9);
	minEnergy.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST_F(ModuleTest, MaximumTrajectoryLength_Continue) {
	MaximumTrajectoryLength maxLength(10);
	candidate.setTrajectoryLength(9.9);
	maxLength.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST_F(ModuleTest, MaximumTrajectoryLength_Break) {
	MaximumTrajectoryLength maxLength(10);
	candidate.setTrajectoryLength(10.1);
	maxLength.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST_F(ModuleTest, SmallObserverSphere) {
	SmallObserverSphere obs(Vector3(0, 0, 0), 1);
	candidate.current.setPosition(Vector3(2, 0, 0));
	obs.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);

	candidate.setStatus(Candidate::Active);
	candidate.current.setPosition(Vector3(0.1, 0.5, -0.1));
	obs.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Detected);

	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(0, 0, 2));
	obs.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1);
}

TEST_F(ModuleTest, CubicBoundary_Inside) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, Candidate::Stopped);
	candidate.current.setPosition(Vector3(9, 5, 5));
	cube.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST_F(ModuleTest, CubicBoundary_Above) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, Candidate::Stopped);
	candidate.current.setPosition(Vector3(10.1, 5, 5));
	cube.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST_F(ModuleTest, CubicBoundary_Below) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, Candidate::Stopped);
	candidate.current.setPosition(Vector3(5, -0.1, 5));
	cube.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST_F(ModuleTest, CubicBoundary_LimitStep) {
	CubicBoundary cube(Vector3(0, 0, 0), 10, 1, Candidate::Stopped);
	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(5, 5, 0.5));
	cube.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST_F(ModuleTest, SphericalBoundary_Inside) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, Candidate::Stopped);
	candidate.current.setPosition(Vector3(9, 0, 0));
	sphere.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST_F(ModuleTest, SphericalBoundary_Outside) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, Candidate::Stopped);
	candidate.current.setPosition(Vector3(0, -10.1, 0));
	sphere.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST_F(ModuleTest, SphericalBoundary_LimitStep) {
	SphericalBoundary sphere(Vector3(0, 0, 0), 10, 1, Candidate::Stopped);
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3(0, 0, 9.5));
	sphere.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

TEST_F(ModuleTest, EllipsoidalBoundary_Inside) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15,
			Candidate::Stopped);
	candidate.current.setPosition(Vector3(3, 2, 0));
	ellipsoid.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
}

TEST_F(ModuleTest, EllipsoidalBoundary_Outside) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15,
			Candidate::Stopped);
	candidate.current.setPosition(Vector3(0, 25, 0));
	ellipsoid.process(&candidate, secondaries);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST_F(ModuleTest, EllipsoidalBoundary_LimitStep) {
	EllipsoidalBoundary ellipsoid(Vector3(-5, 0, 0), Vector3(5, 0, 0), 15, 0.5,
			Candidate::Stopped);
	candidate.setNextStep(2);
	candidate.current.setPosition(Vector3(7, 0, 0));
	ellipsoid.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

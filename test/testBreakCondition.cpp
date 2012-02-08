#include "mpc/module/BreakCondition.h"
#include "mpc/Candidate.h"
#include "gtest/gtest.h"

namespace mpc {

TEST(testBreakCondition, MinimumEnergy) {
	MinimumEnergy minEnergy(5 * EeV);
	Candidate candidate;

	candidate.setStatus(Candidate::Active);
	candidate.current.setEnergy(5.1 * EeV);
	minEnergy.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);

	candidate.setStatus(Candidate::Active);
	candidate.current.setEnergy(4.9 * EeV);
	minEnergy.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testBreakCondition, MaximumTrajectoryLength) {
	MaximumTrajectoryLength maxLength(10 * Mpc);
	Candidate candidate;

	candidate.setStatus(Candidate::Active);
	candidate.setTrajectoryLength(9.9 * Mpc);
	maxLength.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);

	candidate.setStatus(Candidate::Active);
	candidate.setTrajectoryLength(10 * Mpc);
	maxLength.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
}

TEST(testBreakCondition, LargeObserverSphere) {
	LargeObserverSphere obs(Vector3(0, 0, 0), 10);
	Candidate candidate;

	candidate.current.setPosition(Vector3(0, 0, 0));
	obs.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);

	candidate.setStatus(Candidate::Active);
	candidate.current.setPosition(Vector3(10., 0, 0));
	obs.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Detected);

	candidate.setStatus(Candidate::Active);
	candidate.current.setPosition(Vector3(0, -10., 0));
	obs.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Detected);

	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(0, 0, 8));
	obs.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 2);
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

TEST(testBreakCondition, SimulationBox) {
	SimulationBox simBox(Vector3(0, 0, 0), 10, 1);
	Candidate candidate;

	candidate.current.setPosition(Vector3(9, 5, 5));
	simBox.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);

	candidate.setStatus(Candidate::Active);
	candidate.current.setPosition(Vector3(10.1, 5, 5));
	simBox.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::OutOfBounds);

	candidate.setStatus(Candidate::Active);
	candidate.current.setPosition(Vector3(5, -0.1, 1));
	simBox.process(&candidate);
	EXPECT_EQ(candidate.getStatus(), Candidate::OutOfBounds);

	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(9.5, 5, 5));
	simBox.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);

	candidate.setNextStep(10);
	candidate.current.setPosition(Vector3(5, 5, 0.5));
	simBox.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 1.5);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

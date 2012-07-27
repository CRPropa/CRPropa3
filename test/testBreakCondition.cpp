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

TEST(PeriodicBox, high) {
	// Tests if the periodical boundaries place the particle back inside the box and translate the initial position accordingly.
	Vector3d origin(2, 2, 2);
	Vector3d size(2, 2, 2);
	PeriodicBox box(origin, size);

	Candidate c;
	c.current.setPosition(Vector3d(4.5, 4.3, 4.4));
	c.initial.setPosition(Vector3d(3, 3, 3));

	box.process(&c);

	EXPECT_DOUBLE_EQ(2.5, c.current.getPosition().x);
	EXPECT_DOUBLE_EQ(1, c.initial.getPosition().x);
	EXPECT_DOUBLE_EQ(2.3, c.current.getPosition().y);
	EXPECT_DOUBLE_EQ(1, c.initial.getPosition().y);
	EXPECT_DOUBLE_EQ(2.4, c.current.getPosition().z);
	EXPECT_DOUBLE_EQ(1, c.initial.getPosition().z);
}

TEST(PeriodicBox, low) {
	// Tests if the periodical boundaries place the particle back inside the box and translate the initial position accordingly.
	Vector3d origin(0, 0, 0);
	Vector3d size(2, 2, 2);
	PeriodicBox box(origin, size);

	Candidate c;
	c.current.setPosition(Vector3d(-2.5, -0.3, -0.4));
	c.initial.setPosition(Vector3d(1, 1, 1));

	box.process(&c);

	EXPECT_DOUBLE_EQ(1.5, c.current.getPosition().x);
	EXPECT_DOUBLE_EQ(5, c.initial.getPosition().x);
	EXPECT_DOUBLE_EQ(1.7, c.current.getPosition().y);
	EXPECT_DOUBLE_EQ(3, c.initial.getPosition().y);
	EXPECT_DOUBLE_EQ(1.6, c.current.getPosition().z);
	EXPECT_DOUBLE_EQ(3, c.initial.getPosition().z);
}

TEST(ReflectiveBox, high) {
	// Tests if the reflective boundaries place the particle back inside the box and translate the initial position accordingly.
	// Also the initial and final directions are to be reflected
	Vector3d origin(10, 10, 10);
	Vector3d size(10, 20, 20);
	ReflectiveBox box(origin, size);

	Candidate c;
	c.initial.setPosition(Vector3d(15, 15, 15));
	c.initial.setDirection(Vector3d(0, 0.6, 0.8));
	c.current.setPosition(Vector3d(15, 15, 31));
	c.current.setDirection(Vector3d(0, 0.6, 0.8));

	box.process(&c);

	EXPECT_DOUBLE_EQ(15, c.initial.getPosition().x);
	EXPECT_DOUBLE_EQ(15, c.initial.getPosition().y);
	EXPECT_DOUBLE_EQ(45, c.initial.getPosition().z);

	EXPECT_DOUBLE_EQ(15, c.current.getPosition().x);
	EXPECT_DOUBLE_EQ(15, c.current.getPosition().y);
	EXPECT_DOUBLE_EQ(29, c.current.getPosition().z);

	EXPECT_DOUBLE_EQ(0, c.initial.getDirection().x);
	EXPECT_DOUBLE_EQ(0.6, c.initial.getDirection().y);
	EXPECT_DOUBLE_EQ(-0.8, c.initial.getDirection().z);

	EXPECT_DOUBLE_EQ(0, c.current.getDirection().x);
	EXPECT_DOUBLE_EQ(0.6, c.current.getDirection().y);
	EXPECT_DOUBLE_EQ(-0.8, c.current.getDirection().z);
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

TEST(CubicBoundary, limitStepLower) {
	CubicBoundary cube(Vector3d(10, 10, 10), 10);
	cube.setLimitStep(true, 1);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(15, 15, 10.5));
	candidate.setNextStep(100);
	cube.process(&candidate);
	EXPECT_DOUBLE_EQ(1.5, candidate.getNextStep());
}

TEST(CubicBoundary, limitStepUpper) {
	CubicBoundary cube(Vector3d(-10, -10, -10), 10);
	cube.setLimitStep(true, 1);
	Candidate candidate;
	candidate.current.setPosition(Vector3d(-5, -5, -0.5));
	candidate.setNextStep(100);
	cube.process(&candidate);
	EXPECT_DOUBLE_EQ(1.5, candidate.getNextStep());
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
	candidate.setNextStep(100);
	candidate.current.setPosition(Vector3d(0, 0, 9.5));
	sphere.process(&candidate);
	EXPECT_DOUBLE_EQ(1.5, candidate.getNextStep());
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

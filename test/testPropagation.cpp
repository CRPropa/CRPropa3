#include "mpc/Candidate.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(testSimplePropagation, step) {
	double minStep = 20;
	double maxStep = 100;
	SimplePropagation propa(minStep, maxStep);

	ParticleState p;
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	c.setNextStep(10);

	propa.process(&c);

	EXPECT_EQ(minStep, c.getCurrentStep());
	EXPECT_EQ(maxStep, c.getNextStep());
	EXPECT_EQ(Vector3d(0,  0, 0), c.initial.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.initial.getDirection());
	EXPECT_EQ(Vector3d(0,  0, 0), c.previous.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.previous.getDirection());
	EXPECT_EQ(Vector3d(0, 20, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.current.getDirection());
}

TEST(testDeflectionCK, proton) {
	DeflectionCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(0.1 * kpc, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(0.5 * kpc, c.getNextStep());
}

TEST(testDeflectionCK, neutron) {
	DeflectionCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(0.1 * kpc, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(4000 * Mpc, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 0.1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

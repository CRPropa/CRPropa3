#include "mpc/Candidate.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/magneticField/UniformMagneticField.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(testSimplePropagation, noMinStep) {
	ParticleState p;
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	c.setNextStep(10);

	double acceleration = 5;
	double minStep = 0;
	SimplePropagation propa(acceleration, minStep);
	propa.process(&c);

	EXPECT_EQ(10, c.getCurrentStep());
	EXPECT_EQ(50, c.getNextStep());
	EXPECT_EQ(Vector3d(0,10,0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0,1,0), c.current.getDirection());
}

TEST(testSimplePropagation, withMinStep) {
	ParticleState p;
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	c.setNextStep(10);

	double acceleration = 5;
	double minStep = 20;
	SimplePropagation propa(acceleration, minStep);
	propa.process(&c);

	EXPECT_EQ(20, c.getCurrentStep());
	EXPECT_EQ(100, c.getNextStep());
	EXPECT_EQ(Vector3d(0,20,0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0,1,0), c.current.getDirection());
}

TEST(testDeflectionCK, proton) {
	ParticleState p;
	p.setId(getNucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	c.setNextStep(0);

	ref_ptr<UniformMagneticField> field = new UniformMagneticField(
			Vector3d(0, 0, 1e-12));

	DeflectionCK propa(field);
	propa.process(&c);

	EXPECT_DOUBLE_EQ(0.1 * kpc, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(0.5 * kpc, c.getNextStep());
}

TEST(testDeflectionCK, neutron) {
	ParticleState p;
	p.setId(getNucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);

	ref_ptr<UniformMagneticField> field = new UniformMagneticField(
			Vector3d(0, 0, 1e-12));

	DeflectionCK propa(field);
	propa.process(&c);

	EXPECT_DOUBLE_EQ(0.1 * kpc, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(0.5 * kpc, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 0.1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

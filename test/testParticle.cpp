#include "mpc/Particle.h"
#include "gtest/gtest.h"

using mpc::Mpc;
using mpc::EeV;

TEST(particleTest, position) {
	mpc::Particle particle;
	mpc::Hep3Vector v(1, 3, 5);
	particle.setPosition(v * Mpc);

	EXPECT_DOUBLE_EQ(particle.getPosition().getX() / Mpc, 1);
	EXPECT_DOUBLE_EQ(particle.getPosition().getY() / Mpc, 3);
	EXPECT_DOUBLE_EQ(particle.getPosition().getZ() / Mpc, 5);

	EXPECT_DOUBLE_EQ(particle.getPosition().getX(), 1 * Mpc);
	EXPECT_DOUBLE_EQ(particle.getPosition().getY(), 3 * Mpc);
	EXPECT_DOUBLE_EQ(particle.getPosition().getZ(), 5 * Mpc);
}

TEST(particleTest, direction) {
	mpc::Particle particle;
	mpc::Hep3Vector v(1, 2, 3);
	particle.setDirection(v);
	particle.setEnergy(10 * EeV);

	EXPECT_EQ(particle.getEnergy(), 10* EeV);
	EXPECT_EQ(particle.getEnergy(), 10*EeV);

	EXPECT_TRUE(particle.getDirection() == v);
	EXPECT_TRUE(particle.getVelocity() == v * mpc::c_light);
	EXPECT_TRUE(
			particle.getMomentum() == v * (particle.getEnergy() / mpc::c_light));
}

TEST(particleTest, proton) {
	mpc::Particle particle;
	particle.setChargeNumber(1);
	particle.setMass(1 * mpc::amu);

	EXPECT_EQ(particle.getChargeNumber(), 1);
	EXPECT_EQ(particle.getMass(), 1 * mpc::amu);
	EXPECT_DOUBLE_EQ(particle.getChargeNumber(), 1);
//	EXPECT_DOUBLE_EQ(particle.getMass(), 1.660538921e-27);

	particle.setEnergy(mpc::amu * mpc::c_squared);
	EXPECT_DOUBLE_EQ(particle.getLorentzFactor(), 1);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

#include "Particle.h"
#include "gtest/gtest.h"

TEST(particleTest, position) {
	Particle particle;
	Hep3Vector v(1,3,5);
	particle.setPositionMpc(v);

	EXPECT_DOUBLE_EQ(particle.getPositionMpc().getX(), 1);
	EXPECT_DOUBLE_EQ(particle.getPositionMpc().getY(), 3);
	EXPECT_DOUBLE_EQ(particle.getPositionMpc().getZ(), 5);

	EXPECT_DOUBLE_EQ(particle.getPosition().getX(), 1*Mpc);
	EXPECT_DOUBLE_EQ(particle.getPosition().getY(), 3*Mpc);
	EXPECT_DOUBLE_EQ(particle.getPosition().getZ(), 5*Mpc);
}

TEST(particleTest, direction) {
	Particle particle;
	Hep3Vector v(1,2,3);
	particle.setDirection(v);
	particle.setEnergyEeV(10);

	EXPECT_EQ(particle.getEnergyEeV(), 10);
	EXPECT_EQ(particle.getEnergy(), 10*EeV);

	EXPECT_TRUE(particle.getDirection() == v);
	EXPECT_TRUE(particle.getVelocity() == v * c_light);
	EXPECT_TRUE(particle.getMomentum() == v * (particle.getEnergy() / c_light));
}

TEST(particleTest, proton) {
	Particle particle;
	particle.setChargeMassNumber(1,1);

	EXPECT_EQ(particle.getChargeNumber(), 1);
	EXPECT_EQ(particle.getMassNumber(), 1);
	EXPECT_DOUBLE_EQ(particle.getCharge(), 1.602176487e-19);
	EXPECT_DOUBLE_EQ(particle.getMass(), 1.660538921e-27);

	particle.setEnergy(amu * c_squared);
	EXPECT_DOUBLE_EQ(particle.getLorentzFactor(), 1);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#include "mpc/Candidate.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(ParticleState, position) {
	ParticleState particle;
	Vector3 v(1, 3, 5);
	particle.setPosition(v * Mpc);
	EXPECT_TRUE(particle.getPosition() == v * Mpc);
}

TEST(ParticleState, energy) {
	ParticleState particle;
	particle.setEnergy(10 * EeV);
	EXPECT_EQ(particle.getEnergy(), 10 * EeV);
}

TEST(ParticleState, direction) {
	ParticleState particle;
	Vector3 v(1, 2, 3);
	particle.setDirection(v);
	EXPECT_TRUE(particle.getDirection() == v / v.mag());
}

TEST(ParticleState, velocity) {
	ParticleState particle;
	Vector3 v(1, 1, 0);
	particle.setDirection(v);
	EXPECT_TRUE(particle.getVelocity() == v / v.mag() * c_light);
}

TEST(ParticleState, momentum) {
	ParticleState particle;
	Vector3 v(0, 1, 0);
	particle.setDirection(v);
	particle.setEnergy(100 * EeV);
	EXPECT_TRUE(particle.getMomentum() == v * (particle.getEnergy() / c_light));
}

TEST(ParticleState, id) {
	ParticleState particle;
	particle.setId(getNucleusId(12, 6));
	EXPECT_EQ(particle.getId(), 1000060120);
}

TEST(ParticleState, idException) {
	EXPECT_THROW(getNucleusId(5, 6), std::runtime_error);
}

TEST(ParticleState, charge) {
	ParticleState particle;
	particle.setId(getNucleusId(56, 26));
	EXPECT_EQ(particle.getChargeNumber(), 26);
	EXPECT_DOUBLE_EQ(particle.getCharge(), 26 * eplus);
}

TEST(ParticleState, massProton) {
	ParticleState particle;
	particle.setId(getNucleusId(1, 1));
	EXPECT_EQ(particle.getMassNumber(), 1);
	EXPECT_DOUBLE_EQ(particle.getMass(), mass_proton);
}

TEST(ParticleState, massNeutron) {
	ParticleState particle;
	particle.setId(getNucleusId(1, 0));
	EXPECT_EQ(particle.getMassNumber(), 1);
	EXPECT_DOUBLE_EQ(particle.getMass(), mass_neutron);
}

TEST(ParticleState, lorentzFactor) {
	ParticleState particle;
	particle.setId(getNucleusId(1, 1));
	particle.setEnergy(1e12 * eV);
	EXPECT_DOUBLE_EQ(particle.getLorentzFactor(), 1e12 * eV / mass_proton / c_squared);
}

TEST(Candidate, currentStep) {
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getCurrentStep(), 1 * Mpc);
}

TEST(Candidate, limitNextStep) {
	Candidate candidate;
	candidate.setNextStep(5 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 5 * Mpc);
	candidate.limitNextStep(2 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 2 * Mpc);
	candidate.limitNextStep(3 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 2 * Mpc);
}

TEST(Candidate, isActive) {
	Candidate candidate;
	EXPECT_TRUE(candidate.isActive());
	candidate.setActive(false);
	EXPECT_FALSE(candidate.isActive());
}

TEST(Candidate, property) {
	Candidate candidate;
	candidate.setProperty("foo", "bar");
	EXPECT_TRUE(candidate.hasProperty("foo"));
	std::string value;
	candidate.getProperty("foo", value);
	EXPECT_EQ("bar", value);
}

TEST(common, digit) {
	EXPECT_EQ(1, digit(1234, 1000));
	EXPECT_EQ(2, digit(1234, 100));
	EXPECT_EQ(3, digit(1234, 10));
	EXPECT_EQ(4, digit(1234, 1));
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

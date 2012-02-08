#include "mpc/Candidate.h"
#include "gtest/gtest.h"

namespace mpc {

TEST(testParticleState, position) {
	ParticleState particle;
	Vector3 v(1, 3, 5);
	particle.setPosition(v * Mpc);
	EXPECT_TRUE(particle.getPosition() == v * Mpc);
}

TEST(testParticleState, energy) {
	ParticleState particle;
	particle.setEnergy(10 * EeV);
	EXPECT_EQ(particle.getEnergy(), 10 * EeV);
}

TEST(testParticleState, direction) {
	ParticleState particle;
	Vector3 v(1, 2, 3);
	particle.setDirection(v);
	EXPECT_TRUE(particle.getDirection() == v / v.mag());
}

TEST(testParticleState, velocity) {
	ParticleState particle;
	Vector3 v(1, 1, 0);
	particle.setDirection(v);
	EXPECT_TRUE(particle.getVelocity() == v / v.mag() * c_light);
}

TEST(testParticleState, momentum) {
	ParticleState particle;
	Vector3 v(0, 1, 0);
	particle.setDirection(v);
	particle.setEnergy(100 * EeV);
	EXPECT_TRUE(particle.getMomentum() == v * (particle.getEnergy() / c_light));
}

TEST(testParticleState, id) {
	ParticleState particle;
	particle.setId(1045026000);
	EXPECT_EQ(particle.getId(), 1045026000);
}

TEST(testParticleState, charge) {
	ParticleState particle;
	particle.setId(1056026000);
	EXPECT_EQ(particle.getChargeNumber(), 26);
	EXPECT_DOUBLE_EQ(particle.getCharge(), 26 * eplus);
}

TEST(testParticleState, mass) {
	ParticleState particle;
	particle.setId(1056026000);
	EXPECT_EQ(particle.getMassNumber(), 56);
	EXPECT_DOUBLE_EQ(particle.getMass(), 56 * amu);
}

TEST(testParticleState, lorentzFactor) {
	ParticleState particle;
	particle.setId(1010005000);
	particle.setEnergy(1e12 * eV);
	double lf = 1e12 * eV / (10 * amu * c_squared);
	EXPECT_DOUBLE_EQ(particle.getLorentzFactor(), lf);
}

TEST(testCandidate, currentStep) {
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getCurrentStep(), 1 * Mpc);
}

TEST(testCandidate, limitNextStep) {
	Candidate candidate;
	candidate.setNextStep(5 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 5 * Mpc);
	candidate.limitNextStep(2 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 2 * Mpc);
	candidate.limitNextStep(3 * Mpc);
	EXPECT_DOUBLE_EQ(candidate.getNextStep(), 2 * Mpc);
}

TEST(testCandidate, status) {
	Candidate candidate;
	candidate.setStatus(Candidate::Active);
	EXPECT_EQ(candidate.getStatus(), Candidate::Active);
	candidate.setStatus(Candidate::Detected);
	EXPECT_EQ(candidate.getStatus(), Candidate::Detected);
	candidate.setStatus(Candidate::Stopped);
	EXPECT_EQ(candidate.getStatus(), Candidate::Stopped);
	candidate.setStatus(Candidate::UserDefined);
	EXPECT_EQ(candidate.getStatus(), Candidate::UserDefined);
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

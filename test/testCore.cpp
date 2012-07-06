#include "mpc/Candidate.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(ParticleState, position) {
	ParticleState particle;
	Vector3d v(1, 3, 5);
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
	Vector3d v(1, 2, 3);
	particle.setDirection(v);
	EXPECT_TRUE(particle.getDirection() == v.getUnitVector());
}

TEST(ParticleState, velocity) {
	ParticleState particle;
	Vector3d v(1, 1, 0);
	particle.setDirection(v);
	EXPECT_TRUE(particle.getVelocity() == v.getUnitVector() * c_light);
}

TEST(ParticleState, momentum) {
	ParticleState particle;
	Vector3d v(0, 1, 0);
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
	EXPECT_DOUBLE_EQ(particle.getLorentzFactor(),
			1e12 * eV / mass_proton / c_squared);
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

TEST(common, interpolate) {
	double xD[100];
	double yD[100];
	for (int i = 0; i < 100; i++) {
		xD[i] = 1 + i * 0.02;
		yD[i] = pow(xD[i], 2);
	}
	double yInt = interpolate(1.5001, xD, yD);
	EXPECT_NEAR(pow(1.5001, 2), yInt, 1e-4);
}

TEST(common, interpolateVector) {
	std::vector<double> xD, yD;
	xD.resize(100);
	yD.resize(100);
	for (int i = 0; i < 100; i++) {
		xD[i] = 1 + i * 0.02;
		yD[i] = pow(xD[i], 2);
	}
	double yInt = interpolate(1.5001, &xD[0], &yD[0]);
	EXPECT_NEAR(pow(1.5001, 2), yInt, 1e-4);
}

TEST(common, interpolateEquidistant) {
	double yD[100];
	for (int i = 0; i < 100; i++) {
		yD[i] = pow(1 + i * 0.02, 2);
	}
	double yInt = interpolateEquidistant(1.5001, 1, 0.02, yD);
	EXPECT_NEAR(pow(1.5001, 2), yInt, 1e-4);
}

TEST(Vector3, division) {
	Vector3d v(10);
	v /= 10;
	EXPECT_DOUBLE_EQ(v.x, 1);
	EXPECT_DOUBLE_EQ(v.y, 1);
	EXPECT_DOUBLE_EQ(v.z, 1);
	v = Vector3d(10) / Vector3d(2); // element-wise division
	EXPECT_DOUBLE_EQ(v.x, 5);
	EXPECT_DOUBLE_EQ(v.y, 5);
	EXPECT_DOUBLE_EQ(v.z, 5);
}

TEST(Vector3, mod) {
	Vector3d v(10.1, 10.2, 10.3);
	v %= 10.2;
	EXPECT_NEAR(v.x, 10.1, 1e-10); // mod doesn't preserve double precision
	EXPECT_NEAR(v.y, 0, 1e-10);
	EXPECT_NEAR(v.z, 0.1, 1e-10);
}

TEST(Vector3, dot) {
	double a = Vector3d(1, 0, 0).dot(Vector3d(0, 1, 0));
	EXPECT_DOUBLE_EQ(a, 0);
	double b = Vector3d(-1, 10, 2).dot(Vector3d(5, 1, -3));
	EXPECT_DOUBLE_EQ(b, -1);
}

TEST(Vector3, cross) {
	Vector3d v = Vector3d(1, 0, 0).cross(Vector3d(0, 1, 0));
	EXPECT_DOUBLE_EQ(v.x, 0);
	EXPECT_DOUBLE_EQ(v.y, 0);
	EXPECT_DOUBLE_EQ(v.z, 1);
}

TEST(Vector3, angle) {
	double a = Vector3d(1, 1, 0).getAngleTo(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(a, 45 * M_PI / 180);
	double b = Vector3d(0, 0, 1).getAngleTo(Vector3d(0, 0, 1));
	EXPECT_DOUBLE_EQ(b, 0);
}

TEST(Vector3, magnitude) {
	Vector3d v = Vector3d(1, 2, -2);
	EXPECT_DOUBLE_EQ(v.getMag(), 3);
	EXPECT_DOUBLE_EQ(v.getMag2(), 9);
}

TEST(Vector3, distance) {
	double a = Vector3d(10, 0, 10).getDistanceTo(Vector3d(10, 0, 0));
	EXPECT_DOUBLE_EQ(a, 10);
}

TEST(Vector3, rotation) {
	Vector3d v(10, 0, 0);
	v.rotate(Vector3d(0, 2, 0), M_PI);
	EXPECT_NEAR(v.x, -10, 1e-9);
	// rotation doesn't preserve double precision
	EXPECT_NEAR(v.y, 0, 1e-9);
	EXPECT_NEAR(v.z, 0, 1e-9);

	v = Vector3d(10, 0, 0);
	v.rotate(Vector3d(1, 0, 0), M_PI); // no rotation
	EXPECT_NEAR(v.x, 10, 1e-9);
	EXPECT_NEAR(v.y, 0, 1e-9);
	EXPECT_NEAR(v.z, 0, 1e-9);

	v = Vector3d(5, 2, 7);
	v.rotate(Vector3d(1, 8, -4), 2 * M_PI); // full rotation
	EXPECT_NEAR(v.x, 5, 1e-9);
	EXPECT_NEAR(v.y, 2, 1e-9);
	EXPECT_NEAR(v.z, 7, 1e-9);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

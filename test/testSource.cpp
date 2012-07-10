#include "mpc/Source.h"

#include "gtest/gtest.h"

#include <stdexcept>

namespace mpc {

TEST(testSourcePosition, simpleTest) {
	Vector3d position(1, 2, 3);
	SourcePosition source(position);
	ParticleState ps;
	source.prepare(ps);
	EXPECT_EQ(position, ps.getPosition());
}

TEST(testSourceSphericalVolume, simpleTest) {
	Vector3d center(0, 0, 0);
	double radius = 110;
	SourceSphericalVolume source(center, radius);
	ParticleState ps;
	source.prepare(ps);
	double distance = ps.getPosition().getDistanceTo(center);
	EXPECT_GE(radius, distance);
}

TEST(testSourcePowerLawSpectrum, simpleTest) {
	double Emin = 4 * EeV;
	double Emax = 200 * EeV;
	double index = -2.7;
	SourcePowerLawSpectrum spectrum(Emin, Emax, index);
	ParticleState ps;
	spectrum.prepare(ps);
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(testSourceComposition, simpleTest) {
	Vector3d position(1, 2, 3);
	double Emin = 10;
	double Emax = 100;
	double index = -2;
	SourceComposition source(Emin, Emax, index);
	source.add(getNucleusId(6, 3), 1);
	ParticleState ps;
	source.prepare(ps);
	EXPECT_EQ(6, ps.getMassNumber());
	EXPECT_EQ(3, ps.getChargeNumber());
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(testSourceComposition, throwNoIsotope) {
	SourceComposition source(1, 10, -1);
	ParticleState ps;
	EXPECT_THROW(source.prepare(ps), std::runtime_error);
}

TEST(testSource, simpleTest) {
	Source source;
	source.addProperty(new SourcePosition(Vector3d(10, 0, 0) * Mpc));
	source.addProperty(new SourceIsotropicEmission());
	source.addProperty(new SourcePowerLawSpectrum(5 * EeV, 100 * EeV, -2));
	source.addProperty(new SourceParticleType(getNucleusId(8, 4)));
	ParticleState ps;
	source.prepare(ps);
	EXPECT_EQ(8, ps.getMassNumber());
	EXPECT_EQ(4, ps.getChargeNumber());
	EXPECT_LE(5 * EeV, ps.getEnergy());
	EXPECT_GE(100 * EeV , ps.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, ps.getPosition());
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

#include "mpc/Source.h"

#include "gtest/gtest.h"

#include <stdexcept>

namespace mpc {

TEST(testBasicSource, simpleTest) {
	Vector3d position(1, 2, 3);
	int id = getNucleusId(4,2);
	double Emin = 4 * EeV;
	double Emax = 200 * EeV;
	double breakpoint = 25 * EeV;
	double index1 = -2.7;
	double index2 = -4.2;
	BasicSource source(position, id, Emin, Emax, index1, index2, breakpoint);
	ParticleState ps;
	source.prepare(ps);
	EXPECT_EQ(4, ps.getMassNumber());
	EXPECT_EQ(2, ps.getChargeNumber());
	EXPECT_EQ(position, ps.getPosition());
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(testCompositeSource, simpleTest) {
	Vector3d position(1, 2, 3);
	double Emin = 10;
	double Emax = 100;
	double index = -2;
	CompositeSource source(position, Emin, Emax, index);
	source.addToComposition(getNucleusId(6, 3), 1);
	ParticleState ps;
	source.prepare(ps);
	EXPECT_EQ(6, ps.getMassNumber());
	EXPECT_EQ(3, ps.getChargeNumber());
	EXPECT_EQ(position, ps.getPosition());
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(testCompositeSource, spectralIndexOne) {
	Vector3d position(1, 2, 3);
	double Emin = 10;
	double Emax = 100;
	double index = -1;
	CompositeSource source(position, Emin, Emax, index);
	source.addToComposition(getNucleusId(1, 1), 1);
	ParticleState ps;
	source.prepare(ps);
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(testCompositeSource, throwNoIsotope) {
	CompositeSource source;
	ParticleState ps;
	EXPECT_THROW(source.prepare(ps), std::runtime_error);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

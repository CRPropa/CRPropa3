#include "mpc/Source.h"

#include "gtest/gtest.h"
#include <stdexcept>

namespace mpc {

TEST(SourcePosition, simpleTest) {
	Vector3d position(1, 2, 3);
	SourcePosition source(position);
	ParticleState ps;
	source.prepare(ps);
	EXPECT_EQ(position, ps.getPosition());
}

TEST(SourceUniformDistributionSphere, simpleTest) {
	Vector3d center(0, 0, 0);
	double radius = 110;
	SourceUniformDistributionSphere source(center, radius);
	ParticleState ps;
	source.prepare(ps);
	double distance = ps.getPosition().getDistanceTo(center);
	EXPECT_GE(radius, distance);
}

TEST(SourceUniformDistributionBox, simpleTest) {
	Vector3d origin(-7, -2, 0);
	Vector3d size(13, 55, 192);
	SourceUniformDistributionBox box(origin, size);
	ParticleState p;
	box.prepare(p);
	Vector3d pos = p.getPosition();
	EXPECT_LE(origin.x, pos.x);
	EXPECT_LE(origin.y, pos.y);
	EXPECT_LE(origin.z, pos.z);
	EXPECT_GE(size.x, pos.x);
	EXPECT_GE(size.y, pos.y);
	EXPECT_GE(size.z, pos.z);
}

TEST(SourceDensityGrid, withInRange) {
	// Create a grid with 10^3 cells ranging from (0, 0, 0) to (10, 10, 10)
	Vector3d origin(0.5, 0.5, 0.5);
	int cells = 10;
	double spacing = 1;
	ref_ptr<ScalarGrid> grid = new ScalarGrid(origin, cells, spacing);
	for (int ix = 0; ix < cells; ix++)
		for (int iy = 0; iy < cells; iy++)
			for (int iz = 0; iz < cells; iz++)
				grid->get(ix, iy, iz) = ix * iy * iz;

	SourceDensityGrid source(grid);
	ParticleState p;

	source.prepare(p);
	Vector3d pos = p.getPosition();

	// dialed positions should be within the volume (0, 0, 0) - (10, 10, 10)
	EXPECT_LE(0, pos.x);
	EXPECT_GE(10, pos.x);
	EXPECT_LE(0, pos.y);
	EXPECT_GE(10, pos.y);
	EXPECT_LE(0, pos.z);
	EXPECT_GE(10, pos.z);
}

TEST(SourceDensityGrid, OneAllowedCell) {
	// Create a grid with 2^3 cells ranging from (0, 0, 0) to (4, 4, 4)
	Vector3d origin(1, 1, 1);
	int cells = 2;
	double spacing = 2;
	ref_ptr<ScalarGrid> grid = new ScalarGrid(origin, cells, spacing);

	// set all but one cells to 0
	for (int ix = 0; ix < cells; ix++)
		for (int iy = 0; iy < cells; iy++)
			for (int iz = 0; iz < cells; iz++)
				grid->get(ix, iy, iz) = 0;

	// set the first cell ((0, 0, 0) to (2, 2, 2))
	grid->get(0, 0, 0) = 1;

	SourceDensityGrid source(grid);
	ParticleState p;

	int nFalse = 0;
	Vector3d mean(0, 0, 0);
	for (int i = 0; i < 10000; i++) {
		source.prepare(p);
		Vector3d pos = p.getPosition();
		mean += pos;
		if ((pos.x < 0) or (pos.x > 2) or (pos.y < 0) or (pos.y > 2)
				or (pos.z < 0) or (pos.z > 2))
			nFalse++;
	}

	// only the first bin should get dialed
	EXPECT_EQ(0, nFalse);

	// mean should be close to (1, 1, 1) if random positions are uniform in (0, 0, 0) - (2, 2, 2)
	mean /= 10000;
	EXPECT_NEAR(1, mean.x, 0.1);
	EXPECT_NEAR(1, mean.y, 0.1);
	EXPECT_NEAR(1, mean.z, 0.1);
}

TEST(SourceDensityGrid1D, withInRange) {
	// Create a grid with 10 cells ranging from 0 to 10
	Vector3d origin(0.5, 0, 0);
	ref_ptr<ScalarGrid> grid = new ScalarGrid(origin, 10, 1, 1, 1.0);

	for (int i = 0; i < 10; i++) {
		grid->get(i, 0, 0) = 2;
	}

	SourceDensityGrid1D source(grid);
	ParticleState p;

	source.prepare(p);
	Vector3d pos = p.getPosition();
	// dialed position should be within the range 0 - 10
	EXPECT_LE(0, pos.x);
	EXPECT_GE(10, pos.x);
}

TEST(SourceDensityGrid1D, OneAllowedCell) {
	Vector3d origin(0.5, 0, 0);
	ref_ptr<ScalarGrid> grid = new ScalarGrid(origin, 10, 1, 1, 1.0);
	for (int i = 0; i < 10; i++) {
		grid->get(i, 0, 0) = 2;
	}

	grid->get(5, 0, 0);

	SourceDensityGrid1D source(grid);
	ParticleState p;

	source.prepare(p);
	Vector3d pos = p.getPosition();
	// dialed position should be within the range 0 - 10
	EXPECT_LE(0, pos.x);
	EXPECT_GE(10, pos.x);
}

TEST(SourcePowerLawSpectrum, simpleTest) {
	double Emin = 4 * EeV;
	double Emax = 200 * EeV;
	double index = -2.7;
	SourcePowerLawSpectrum spectrum(Emin, Emax, index);
	ParticleState ps;
	spectrum.prepare(ps);

	// energy should be within Emin - Emax
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(SourceComposition, simpleTest) {
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

TEST(SourceComposition, throwNoIsotope) {
	SourceComposition source(1, 10, -1);
	ParticleState ps;
	EXPECT_THROW(source.prepare(ps), std::runtime_error);
}

TEST(Source, allPropertiesUsed) {
	Source source;
	source.addProperty(new SourcePosition(Vector3d(10, 0, 0) * Mpc));
	source.addProperty(new SourceIsotropicEmission());
	source.addProperty(new SourcePowerLawSpectrum(5 * EeV, 100 * EeV, -2));
	source.addProperty(new SourceParticleType(getNucleusId(8, 4)));

	Candidate c = *source.getCandidate();

	ParticleState p = c.initial;
	EXPECT_EQ(8, p.getMassNumber());
	EXPECT_EQ(4, p.getChargeNumber());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());

	p = c.previous;
	EXPECT_EQ(8, p.getMassNumber());
	EXPECT_EQ(4, p.getChargeNumber());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());

	p = c.current;
	EXPECT_EQ(8, p.getMassNumber());
	EXPECT_EQ(4, p.getChargeNumber());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());
}

TEST(SourceList, simpleTest) {
	// test if source list works with one source
	SourceList sourceList;
	ref_ptr<Source> source = new Source;
	source->addProperty(new SourcePosition(Vector3d(10, 0, 0)));
	sourceList.addSource(source);

	ref_ptr<Candidate> c = sourceList.getCandidate();

	EXPECT_EQ(Vector3d(10, 0, 0), c->initial.getPosition());
	EXPECT_EQ(Vector3d(10, 0, 0), c->previous.getPosition());
	EXPECT_EQ(Vector3d(10, 0, 0), c->current.getPosition());
}

TEST(SourceList, noSource) {
	// test if an error is thrown when source list empty
	SourceList sourceList;
	EXPECT_THROW(sourceList.getCandidate(), std::runtime_error);
}

TEST(SourceList, luminosity) {
	// test if the sources are dialed according to their luminosities
	SourceList sourceList;

	ref_ptr<Source> source1 = new Source;
	source1->addProperty(new SourceEnergy(100));
	sourceList.addSource(source1, 80);

	ref_ptr<Source> source2 = new Source;
	source2->addProperty(new SourceEnergy(0));
	sourceList.addSource(source2, 20);

	double meanE = 0;
	for (int i = 0; i < 1000; i++) {
		ref_ptr<Candidate> c = sourceList.getCandidate();
		meanE += c->initial.getEnergy();
	}
	meanE /= 1000;
	EXPECT_NEAR(80, meanE, 2); // this test can stochastically fail
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

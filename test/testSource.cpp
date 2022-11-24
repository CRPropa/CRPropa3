#include "crpropa/Source.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"

#include "gtest/gtest.h"
#include <stdexcept>

namespace crpropa {

TEST(SourcePosition, simpleTest) {
	Vector3d position(1, 2, 3);
	SourcePosition source(position);
	ParticleState ps;
	source.prepareParticle(ps);
	EXPECT_EQ(position, ps.getPosition());
}

TEST(SourceMultiplePositions, simpleTest) {
	SourceMultiplePositions source;
	source.add(Vector3d(1, 0, 0), 0.25);
	source.add(Vector3d(2, 0, 0), 0.75);
	ParticleState ps;
	int n1 = 0;
	int n2 = 0;
	for (int i = 0; i < 10000; i++) {
		source.prepareParticle(ps);
		if (ps.getPosition().x == 1)
			n1++;
		else if (ps.getPosition().x == 2)
			n2++;
	}
	EXPECT_NEAR(n1, 2500, 5 * sqrt(2500));
	EXPECT_NEAR(n2, 7500, 5 * sqrt(7500));
}

TEST(SourceUniformSphere, simpleTest) {
	Vector3d center(0, 0, 0);
	double radius = 110;
	SourceUniformSphere source(center, radius);
	ParticleState ps;
	source.prepareParticle(ps);
	double distance = ps.getPosition().getDistanceTo(center);
	EXPECT_GE(radius, distance);
}

TEST(SourceUniformHollowSphere, simpleTest) {
	Vector3d center(0, 0, 0);
	double radius_inner = 50;
	double radius_outer = 110;
	SourceUniformHollowSphere source(center,
			radius_inner,
			radius_outer);
	for (int i=0; i < 100; ++i) {
		ParticleState ps;
		source.prepareParticle(ps);
		double distance = ps.getPosition().getDistanceTo(center);
		EXPECT_GE(radius_outer, distance);
		EXPECT_LE(radius_inner, distance);
	}
}

TEST(SourceUniformBox, simpleTest) {
	Vector3d origin(-7, -2, 0);
	Vector3d size(13, 55, 192);
	SourceUniformBox box(origin, size);
	ParticleState p;
	box.prepareParticle(p);
	Vector3d pos = p.getPosition();
	EXPECT_LE(origin.x, pos.x);
	EXPECT_LE(origin.y, pos.y);
	EXPECT_LE(origin.z, pos.z);
	EXPECT_GE(size.x, pos.x);
	EXPECT_GE(size.y, pos.y);
	EXPECT_GE(size.z, pos.z);
}

TEST(SourceUniformCylinder, simpleTest) {
	Vector3d center(0, 0, 0);
	double radius = 15;
	double height = 2;
	SourceUniformCylinder cylinder(center, height, radius);
	ParticleState ps;
	cylinder.prepareParticle(ps);
	Vector3d pos = ps.getPosition();
	double R2 = pos.x*pos.x+pos.y*pos.y;
	double H = pow(pos.z*pos.z, 0.5);
	EXPECT_GE(radius*radius, R2);
	EXPECT_GE(height/2., H);
}

TEST(SourceSNRDistribution, simpleTest) {
	double R_earth = 8.5*kpc;
	double alpha = 2.0;
	double beta = 3.53;
	double Z_G = 0.3*kpc;
	SourceSNRDistribution snr(R_earth,alpha, beta, Z_G);
	ParticleState ps;
	snr.prepareParticle(ps);
	Vector3d pos = ps.getPosition();
	double R2 = pos.x*pos.x+pos.y*pos.y;
	EXPECT_GE(20*kpc*20*kpc, R2); // radius must be smaller than 20 kpc
	
	double R2_mean = 0.;
	double Z_mean = 0.;
	for (size_t i=0; i<100000; i++) {
		snr.prepareParticle(ps);
		Vector3d pos = ps.getPosition();
		R2_mean += pow(pos.x/kpc, 2.)+pow(pos.y/kpc, 2.);
		Z_mean += pos.z/kpc;
	}
	R2_mean/=100000.;
	Z_mean/=100000.;
	EXPECT_NEAR(64.4, R2_mean, 1.);
	EXPECT_NEAR(0., Z_mean, 0.1);
}

TEST(SourceDensityGrid, withInRange) {
	// Create a grid with 10^3 cells ranging from (0, 0, 0) to (10, 10, 10)
	Vector3d origin(0, 0, 0);
	int cells = 10;
	double spacing = 1;
	auto grid = new Grid1f(origin, cells, spacing);
	for (int ix = 0; ix < cells; ix++)
		for (int iy = 0; iy < cells; iy++)
			for (int iz = 0; iz < cells; iz++)
				grid->get(ix, iy, iz) = ix * iy * iz;

	SourceDensityGrid source(grid);
	ParticleState p;

	source.prepareParticle(p);
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
	Vector3d origin(0, 0, 0);
	int cells = 2;
	double spacing = 2;
	auto grid = new Grid1f(origin, cells, spacing);
	
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
		source.prepareParticle(p);
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
	EXPECT_NEAR(1, mean.x, 0.2);
	EXPECT_NEAR(1, mean.y, 0.2);
	EXPECT_NEAR(1, mean.z, 0.2);
}

TEST(SourceDensityGrid1D, withInRange) {
	// Create a grid with 10 cells ranging from 0 to 10
	Vector3d origin(0, 0, 0);
	int nCells = 10;
	double spacing = 1.;
	auto grid = new Grid1f(origin, nCells, 1, 1, spacing);

	// set some values
	for (int i = 0; i < 10; i++) {
		grid->get(i, 0, 0) = 2;
	}

	SourceDensityGrid1D source(grid);
	ParticleState p;

	source.prepareParticle(p);
	Vector3d pos = p.getPosition();
	// dialed position should be within the range 0 - 10
	EXPECT_LE(0, pos.x);
	EXPECT_GE(10, pos.x);
}

TEST(SourceDensityGrid1D, OneAllowedCell) {
	// Test if the only allowed cells is repeatedly selected
	Vector3d origin(0, 0, 0);
	int nCells = 10;
	double spacing = 1.;
	auto grid = new Grid1f(origin, nCells, 1, 1, spacing);

	// set some values
	for (int i = 0; i < 10; i++) {
		grid->get(i, 0, 0) = 0;
	}
	grid->get(5, 0, 0) = 1;

	SourceDensityGrid1D source(grid);
	ParticleState p;

	for (int i = 0; i < 100; i++) {
		source.prepareParticle(p);
		// dialed position should be in range 5-6
		Vector3d pos = p.getPosition();
		EXPECT_LE(5, pos.x);
		EXPECT_GE(6, pos.x);
	}
}

TEST(SourcePowerLawSpectrum, simpleTest) {
	double Emin = 4 * EeV;
	double Emax = 200 * EeV;
	double index = -2.7;
	SourcePowerLawSpectrum spectrum(Emin, Emax, index);
	ParticleState ps;
	spectrum.prepareParticle(ps);

	// energy should be within Emin - Emax
	EXPECT_LE(Emin, ps.getEnergy());
	EXPECT_GE(Emax, ps.getEnergy());
}

TEST(SourceComposition, simpleTest) {
	double Emin = 10;
	double Rmax = 100;
	double index = -1;
	SourceComposition source(Emin, Rmax, index);
	source.add(nucleusId(6, 3), 1);
	ParticleState p;
	source.prepareParticle(p);
	EXPECT_EQ(nucleusId(6, 3), p.getId());
	EXPECT_LE(Emin, p.getEnergy());
	EXPECT_GE(6 * Rmax, p.getEnergy());
}

TEST(SourceDirectedEmission, simpleTest) {
	Vector3d mu(1., 0., 0.);
	double kappa = 1000.;
	SourceDirectedEmission source(mu, kappa);
	Candidate c;
	Vector3d meanDir(0., 0., 0.);
	for (size_t i = 0; i < 1000; i++) {
		source.prepareCandidate(c);
		meanDir += c.source.getDirection();
		double w = c.getWeight();
		EXPECT_GE(w, 0.);
	}
	meanDir /= 1000.;
	EXPECT_NEAR(meanDir.x, 1., 0.1);
	EXPECT_NEAR(meanDir.y, 0., 0.01);
	EXPECT_NEAR(meanDir.z, 0., 0.01);
}

TEST(SourceEmissionCone, simpleTest) {
	Vector3d direction(42., 0., 0.);
	double aperture = 1/42.;
	
	SourceEmissionCone source(direction, aperture);
	
	ParticleState p;
	source.prepareParticle(p);
	double angle = direction.getAngleTo(p.getDirection());
	EXPECT_LE(angle, aperture);
}

#ifdef CRPROPA_HAVE_MUPARSER
TEST(SourceGenericComposition, simpleTest) {
	double Emin = 10;
	double Emax = 100;
	SourceGenericComposition source(Emin, Emax, "E^-2");
	int id1 = nucleusId(6, 3);
	int id2 = nucleusId(12, 6);
	source.add(id1, 1);
	source.add(id2, 10);
	ParticleState p;
	size_t id1Count = 0, id2Count = 0;
	size_t ElowCount = 0, EhighCount = 0;
	size_t n = 100000;
	for (size_t i = 0; i < n; i++) {
		source.prepareParticle(p);
		if (p.getId() == id1)
			id1Count++;
		if (p.getId() == id2)
			id2Count++;
		double e = p.getEnergy();
		if ( (e >= Emin) && (e < 20))
			ElowCount++;
		if ( (e >= 20) && (e <= Emax))
			EhighCount++;

	}
	EXPECT_EQ(n, id1Count + id2Count);
	EXPECT_EQ(n, ElowCount + EhighCount);
	EXPECT_NEAR((float)id1Count/(float)id2Count, 0.1, 0.01);
	EXPECT_NEAR((float)ElowCount/(float)EhighCount, 1.25, 0.1);
}
#endif

TEST(SourceComposition, throwNoIsotope) {
	SourceComposition source(1, 10, -1);
	ParticleState ps;
	EXPECT_THROW(source.prepareParticle(ps), std::runtime_error);
}

TEST(SourceRedshiftEvolution, testInRange) {
	Candidate c;

	double zmin = 0.5;
	double zmax = 2.5;

	// general case: m
	SourceRedshiftEvolution source1(3.2, zmin, zmax);
	for (int i = 0; i < 100; i++) {
		source1.prepareCandidate(c);
		EXPECT_LE(zmin, c.getRedshift());
		EXPECT_GE(zmax, c.getRedshift());
	}

	// general case: m = -1
	SourceRedshiftEvolution source2(-1, zmin, zmax);
	for (int i = 0; i < 100; i++) {
		source2.prepareCandidate(c);
		EXPECT_LE(zmin, c.getRedshift());
		EXPECT_GE(zmax, c.getRedshift());
	}
}

TEST(Source, allPropertiesUsed) {
	Source source;
	source.add(new SourcePosition(Vector3d(10, 0, 0) * Mpc));
	source.add(new SourceIsotropicEmission());
	source.add(new SourcePowerLawSpectrum(5 * EeV, 100 * EeV, -2));
	source.add(new SourceParticleType(nucleusId(8, 4)));
	source.add(new SourceRedshift(2));

	Candidate c = *source.getCandidate();

	EXPECT_EQ(2, c.getRedshift());

	ParticleState p;

	p = c.source;
	EXPECT_EQ(nucleusId(8, 4), p.getId());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());

	p = c.created;
	EXPECT_EQ(nucleusId(8, 4), p.getId());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());

	p = c.previous;
	EXPECT_EQ(nucleusId(8, 4), p.getId());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());

	p = c.current;
	EXPECT_EQ(nucleusId(8, 4), p.getId());
	EXPECT_LE(5 * EeV, p.getEnergy());
	EXPECT_GE(100 * EeV, p.getEnergy());
	EXPECT_EQ(Vector3d(10, 0, 0) * Mpc, p.getPosition());
}

TEST(SourceList, simpleTest) {
	// test if source list works with one source
	SourceList sourceList;
	ref_ptr<Source> source = new Source;
	source->add(new SourcePosition(Vector3d(10, 0, 0)));
	sourceList.add(source);

	ref_ptr<Candidate> c = sourceList.getCandidate();

	EXPECT_EQ(Vector3d(10, 0, 0), c->created.getPosition());
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
	source1->add(new SourceEnergy(100));
	sourceList.add(source1, 80);

	ref_ptr<Source> source2 = new Source;
	source2->add(new SourceEnergy(0));
	sourceList.add(source2, 20);

	double meanE = 0;
	for (int i = 0; i < 1000; i++) {
		ref_ptr<Candidate> c = sourceList.getCandidate();
		meanE += c->created.getEnergy();
	}
	meanE /= 1000;
	EXPECT_NEAR(80, meanE, 4); // this test can stochastically fail
}

TEST(SourceTag, sourceTag) {
	SourceTag tag("mySourceTag");
	Candidate c;
	tag.prepareCandidate(c);
	EXPECT_TRUE(c.getTagOrigin() == "mySourceTag");
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

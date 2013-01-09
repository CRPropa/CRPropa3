/** Unit tests for core features of MPC
  	Candidate
  	ParticleState
  	Random
  	Common functions
 */

#include "mpc/Candidate.h"
#include "mpc/Random.h"
#include "mpc/Grid.h"
#include "mpc/GridTools.h"

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
	particle.setId(nucleusId(12, 6));
	EXPECT_EQ(particle.getId(), 1000060120);
}

TEST(ParticleState, idException) {
	EXPECT_THROW(nucleusId(5, 6), std::runtime_error);
}

TEST(ParticleState, Charge) {
	ParticleState particle;

	particle.setId(nucleusId(56, 26)); // iron
	EXPECT_EQ(26, particle.getChargeNumber());
	EXPECT_DOUBLE_EQ(26 * eplus, particle.getCharge());

	particle.setId(-nucleusId(56, 26)); // anti-iron
	EXPECT_EQ(26, particle.getChargeNumber());
	EXPECT_DOUBLE_EQ(-26 * eplus, particle.getCharge());

	particle.setId(11); // electron
	EXPECT_DOUBLE_EQ(-1 * eplus, particle.getCharge());

	particle.setId(-11); // positron
	EXPECT_DOUBLE_EQ(1 * eplus, particle.getCharge());

	particle.setId(12); // electron neutrino
	EXPECT_DOUBLE_EQ(0, particle.getCharge());

	particle.setId(-12); // electron anti-neutrino
	EXPECT_DOUBLE_EQ(0, particle.getCharge());
}

TEST(ParticleState, Mass) {
	ParticleState particle;

	particle.setId(nucleusId(1, 1)); // proton
	EXPECT_EQ(1, particle.getMassNumber());
	EXPECT_DOUBLE_EQ(mass_proton, particle.getMass());

	particle.setId(nucleusId(1, 0)); // neutron
	EXPECT_EQ(1, particle.getMassNumber());
	EXPECT_DOUBLE_EQ(mass_neutron, particle.getMass());

	int id = nucleusId(56, 26);
	particle.setId(id); // iron
	EXPECT_EQ(56, particle.getMassNumber());
	EXPECT_DOUBLE_EQ(nucleusMass(id), particle.getMass());

	particle.setId(-id); // anti-iron
	EXPECT_EQ(56, particle.getMassNumber());
	EXPECT_DOUBLE_EQ(nucleusMass(-id), particle.getMass());
}

TEST(ParticleState, lorentzFactor) {
	ParticleState particle;
	particle.setId(nucleusId(1, 1));
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

TEST(Candidate, addSecondary) {
	Candidate c;
	c.setRedshift(5);
	c.setTrajectoryLength(23);
	c.initial.setId(nucleusId(56,26));
	c.initial.setEnergy(1000);
	c.initial.setPosition(Vector3d(1,2,3));
	c.initial.setDirection(Vector3d(0,0,1));

	c.addSecondary(nucleusId(1,1), 200);
	Candidate s = *c.secondaries[0];

	EXPECT_EQ(nucleusId(1,1), s.current.getId());
	EXPECT_EQ(200, s.current.getEnergy());

	EXPECT_EQ(5, s.getRedshift());
	EXPECT_EQ(23, s.getTrajectoryLength());
	EXPECT_EQ(1000, s.initial.getEnergy());
	EXPECT_TRUE(Vector3d(1,2,3) == s.initial.getPosition());
	EXPECT_TRUE(Vector3d(0,0,1) == s.initial.getDirection());
}

TEST(common, digit) {
	EXPECT_EQ(1, digit(1234, 1000));
	EXPECT_EQ(2, digit(1234, 100));
	EXPECT_EQ(3, digit(1234, 10));
	EXPECT_EQ(4, digit(1234, 1));
}

TEST(common, interpolate) {
	std::vector<double> xD(100), yD(100);
	for (int i = 0; i < 100; i++) {
		xD[i] = 1 + i * 2. / 99.;
		yD[i] = pow(xD[i], 2);
	}

	// interpolated value should be close to computed
	double y = interpolate(1.5001, xD, yD);
	EXPECT_NEAR(pow(1.5001, 2), y, 1e-4);

	// value out of range, return lower bound
	EXPECT_EQ(1, interpolate(0.9, xD, yD));

	// value out of range, return lower bound
	EXPECT_EQ(9, interpolate(3.1, xD, yD));
}

TEST(common, interpolateEquidistant) {
	std::vector<double> yD(100);
	for (int i = 0; i < 100; i++) {
		yD[i] = pow(1 + i * 2. / 99., 2);
	}

	// interpolated value should be close to computed
	double y = interpolateEquidistant(1.5001, 1, 3, yD);
	EXPECT_NEAR(pow(1.5001, 2), y, 1e-4);

	// value out of range, return lower bound
	EXPECT_EQ(1, interpolateEquidistant(0.9, 1, 3, yD));

	// value out of range, return lower bound
	EXPECT_EQ(9, interpolateEquidistant(3.1, 1, 3, yD));
}

TEST(NucleusId, crpropaScheme) {
	// test conversion to and from the CRPropa2 naming scheme
	EXPECT_EQ(nucleusId(56, 26), convertFromCRPropaId(26056));
	EXPECT_EQ(26056, convertToCRPropaId(nucleusId(56, 26)));
}

TEST(Random, seed) {
	Random &a = Random::instance();
	Random &b = Random::instance();

	a.seed(42);
	double r1 = a.rand();

	a.seed(42);
	double r2 = a.rand();

	a.seed(42);
	double r3 = b.rand();

	// seeding should give same random numbers
	EXPECT_EQ(r1, r2);

	// seeding should work for all instances
	EXPECT_EQ(r1, r3);
}

TEST(Grid, PeriodicClamp) {
	// Test correct determination of lower and upper neighbor
	int lo, hi;

	periodicClamp(23.12, 8, lo, hi);
	EXPECT_EQ(7, lo);
	EXPECT_EQ(0, hi);

	periodicClamp(-23.12, 8, lo, hi);
	EXPECT_EQ(0, lo);
	EXPECT_EQ(1, hi);
}

TEST(ScalarGrid, SimpleTest) {
	// Test construction and parameters
	size_t Nx = 5;
	size_t Ny = 8;
	size_t Nz = 10;
	double spacing = 2.0;
	Vector3d origin(1., 2., 3.);

	ScalarGrid grid(origin, Nx, Ny, Nz, spacing);

	EXPECT_TRUE(origin == grid.getOrigin());
	EXPECT_EQ(Nx, grid.getNx());
	EXPECT_EQ(Ny, grid.getNy());
	EXPECT_EQ(Nz, grid.getNz());
	EXPECT_DOUBLE_EQ(spacing, grid.getSpacing());
	EXPECT_EQ(5 * 8 * 10, grid.getGrid().size());

	// Test index handling: get position of grid point (2, 3, 4)
	size_t some_index = 2 * Ny * Nz + 3 * Nz + 4;
	Vector3d position = origin + Vector3d(2, 3, 4) * spacing;
	EXPECT_EQ(position, grid.positionFromIndex(some_index));

	grid.get(2, 3, 4) = 7;
	EXPECT_FLOAT_EQ(7., grid.getGrid()[some_index]);
	EXPECT_FLOAT_EQ(7., grid.interpolate(position));
}

TEST(ScalarGrid, ClosestValue) {
	ScalarGrid grid(Vector3d(0.), 2, 1);
	grid.get(0, 0, 0) = 1;
	grid.get(0, 0, 1) = 2;
	grid.get(0, 1, 0) = 3;
	grid.get(0, 1, 1) = 4;
	grid.get(1, 0, 0) = 5;
	grid.get(1, 0, 1) = 6;
	grid.get(1, 1, 0) = 7;
	grid.get(1, 1, 1) = 8;

	// Closest value
	EXPECT_FLOAT_EQ(1, grid.closestValue(Vector3d(-0.2,  0, 0.4)));
	EXPECT_FLOAT_EQ(2, grid.closestValue(Vector3d(0.2, 0.1, 0.9)));
	EXPECT_FLOAT_EQ(3, grid.closestValue(Vector3d(0.3, 1.2, 0.2)));
	EXPECT_FLOAT_EQ(7, grid.closestValue(Vector3d(0.6, 0.7, 0.4)));
}

TEST(VectorGrid, Interpolation) {
	// Explicitly test trilinear interpolation
	double spacing = 2.793;
	VectorGrid grid(Vector3d(0.), 3, spacing);
	grid.get(0, 0, 1) = Vector3f(1.7, 0., 0.); // set one value

	Vector3d b;
	b = grid.interpolate(Vector3d(0, 0, 1) * spacing);
	EXPECT_FLOAT_EQ(1.7, b.x);

	b = grid.interpolate(Vector3d(0, 0, 0.9) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9, b.x);

	b = grid.interpolate(Vector3d(0, 0, 1.1) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9, b.x);

	b = grid.interpolate(Vector3d(0, 0.15, 0.9) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.85, b.x);

	b = grid.interpolate(Vector3d(0, 2.15, 0.9) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.15, b.x);
}

TEST(VectordGrid, Scale) {
	// Test scaling a field
	ref_ptr<VectorGrid> grid = new VectorGrid(Vector3d(0.), 3, 1);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid->get(ix, iy, iz) = Vector3f(1, 0, 0);

	scaleGrid(grid, 5);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				EXPECT_FLOAT_EQ(5, grid->interpolate(Vector3d(0.7, 0, 0.1)).x);
}

TEST(VectorGrid, Periodicity) {
	// Test for periodic boundaries: grid(x+a*n) = grid(x)
	size_t n = 3;
	double spacing = 3;
	double size = n * spacing;
	VectorGrid grid(Vector3d(0.), n, spacing);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid.get(ix, iy, iz) = Vector3f(iz + ix, iy * iz, ix - iz * iy);

	Vector3d pos(1.2, 2.3, 0.7);
	Vector3f b = grid.interpolate(pos);
	Vector3f b2 = grid.interpolate(pos + Vector3d(1, 0, 0) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	b2 = grid.interpolate(pos + Vector3d(0, 5, 0) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);

	b2 = grid.interpolate(pos + Vector3d(0, 0, -2) * size);
	EXPECT_FLOAT_EQ(b.x, b2.x);
	EXPECT_FLOAT_EQ(b.y, b2.y);
	EXPECT_FLOAT_EQ(b.z, b2.z);
}

TEST(VectorGrid, DumpLoad) {
	// Dump and load a field grid
	ref_ptr<VectorGrid> grid1 = new VectorGrid(Vector3d(0.), 3, 1);
	ref_ptr<VectorGrid> grid2 = new VectorGrid(Vector3d(0.), 3, 1);

	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid1->get(ix, iy, iz) = Vector3f(1, 2, 3);

	dumpGrid(grid1, "testDump.raw");
	loadGrid(grid2, "testDump.raw");

	for (int ix = 0; ix < 3; ix++) {
		for (int iy = 0; iy < 3; iy++) {
			for (int iz = 0; iz < 3; iz++) {
				Vector3f b1 = grid1->get(ix, iy, iz);
				Vector3f b2 = grid2->get(ix, iy, iz);
				EXPECT_FLOAT_EQ(b1.x, b2.x);
				EXPECT_FLOAT_EQ(b1.y, b2.y);
				EXPECT_FLOAT_EQ(b1.z, b2.z);
			}
		}
	}
}

TEST(VectorGrid, DumpLoadTxt) {
	// Dump and load a field grid
	ref_ptr<VectorGrid> grid1 = new VectorGrid(Vector3d(0.), 3, 1);
	ref_ptr<VectorGrid> grid2 = new VectorGrid(Vector3d(0.), 3, 1);

	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid1->get(ix, iy, iz) = Vector3f(ix, iy, iz);

	dumpGridToTxt(grid1, "testDump.txt", 1e4);
	loadGridFromTxt(grid2, "testDump.txt", 1e-4);

	for (int ix = 0; ix < 3; ix++) {
		for (int iy = 0; iy < 3; iy++) {
			for (int iz = 0; iz < 3; iz++) {
				Vector3f b1 = grid1->get(ix, iy, iz);
				Vector3f b2 = grid2->get(ix, iy, iz);
				EXPECT_FLOAT_EQ(b1.x, b2.x);
				EXPECT_FLOAT_EQ(b1.y, b2.y);
				EXPECT_FLOAT_EQ(b1.z, b2.z);
			}
		}
	}
}

TEST(VectorGrid, Speed) {
	// Dump and load a field grid
	VectorGrid grid(Vector3d(0.), 3, 3);
	for (int ix = 0; ix < 3; ix++)
		for (int iy = 0; iy < 3; iy++)
			for (int iz = 0; iz < 3; iz++)
				grid.get(ix, iy, iz) = Vector3f(1, 2, 3);

	Vector3d b;
	for (int i = 0; i < 100000; i++)
		b = grid.interpolate(Vector3d(i));
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

/** Unit tests for core features of CRPropa
  	Candidate
  	ParticleState
  	Random
  	Common functions
 */

#include "crpropa/Candidate.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"
#include "crpropa/Grid.h"
#include "crpropa/GridTools.h"
#include "crpropa/Geometry.h"
#include "crpropa/EmissionMap.h"

#include <HepPID/ParticleIDMethods.hh>
#include "gtest/gtest.h"

namespace crpropa {

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
	EXPECT_DOUBLE_EQ(26 * eplus, particle.getCharge());

	particle.setId(-nucleusId(56, 26)); // anti-iron
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

TEST(ParticleState, Rigidity) {
	ParticleState particle;

	particle.setId(nucleusId(1, 1)); // proton
	particle.setEnergy(1 * EeV);
	EXPECT_EQ(particle.getRigidity(), 1e18);
}

TEST(ParticleState, Mass) {
	ParticleState particle;

	particle.setId(nucleusId(1, 1)); // proton
	EXPECT_DOUBLE_EQ(mass_proton, particle.getMass());

	particle.setId(nucleusId(1, 0)); // neutron
	EXPECT_DOUBLE_EQ(mass_neutron, particle.getMass());

	int id = nucleusId(56, 26);
	particle.setId(id); // iron
	EXPECT_DOUBLE_EQ(nuclearMass(id), particle.getMass());

	particle.setId(-id); // anti-iron
	EXPECT_DOUBLE_EQ(nuclearMass(-id), particle.getMass());

	// approximation for unkown nucleus A * amu - Z * mass_electron
	int A = 238; int Z = 92; // Uranium92
	EXPECT_DOUBLE_EQ(nuclearMass(A, Z), A*amu - Z*mass_electron);
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
	std::string value = candidate.getProperty("foo");
	EXPECT_EQ("bar", value);
}

TEST(Candidate, addSecondary) {
	Candidate c;
	c.setRedshift(5);
	c.setTrajectoryLength(23);
	c.current.setId(nucleusId(56,26));
	c.current.setEnergy(1000);
	c.current.setPosition(Vector3d(1,2,3));
	c.current.setDirection(Vector3d(0,0,1));

	c.addSecondary(nucleusId(1,1), 200);
	Candidate s = *c.secondaries[0];

	EXPECT_EQ(nucleusId(1,1), s.current.getId());
	EXPECT_EQ(200, s.current.getEnergy());

	EXPECT_EQ(5, s.getRedshift());
	EXPECT_EQ(23, s.getTrajectoryLength());
	EXPECT_EQ(1000, s.created.getEnergy());
	EXPECT_TRUE(Vector3d(1,2,3) == s.created.getPosition());
	EXPECT_TRUE(Vector3d(0,0,1) == s.created.getDirection());
}


TEST(Candidate, serialNumber) {
	Candidate::setNextSerialNumber(42);
	Candidate c;
	EXPECT_EQ(43, c.getSourceSerialNumber());
}

TEST(common, digit) {
	EXPECT_EQ(1, digit(1234, 1000));
	EXPECT_EQ(2, digit(1234, 100));
	EXPECT_EQ(3, digit(1234, 10));
	EXPECT_EQ(4, digit(1234, 1));
}

TEST(common, interpolate) {
	// create vectors x = (0, 0.02, ... 2) and y = 2x + 3 = (3, ... 7)
	std::vector<double> xD(101), yD(101);
	for (int i = 0; i <= 100; i++) {
		xD[i] = i * 0.02;
		yD[i] = 2 * xD[i] + 3;
	}

	// interpolating tabulated values of a linear function should produce exact results
	Random &R = Random::instance();
	double x, ytrue, yinterp;
	for (int i = 0; i < 10000; i++) {
		x = R.rand() * 2; // random value between 0 and 2
		ytrue = 2 * x + 3;
		yinterp = interpolate(x, xD, yD);
		EXPECT_DOUBLE_EQ(ytrue, yinterp);
	}

	// test interpolation in first bin
	x = 0.01;
	ytrue = 2 * x + 3;
	yinterp = interpolate(x, xD, yD);
	EXPECT_DOUBLE_EQ(ytrue, yinterp);

	// test interpolation in last bin
	x = 1.99;
	ytrue = 2 * x + 3;
	yinterp = interpolate(x, xD, yD);
	EXPECT_DOUBLE_EQ(ytrue, yinterp);

	// value out of range, return lower bound
	EXPECT_EQ(3, interpolate(-0.001, xD, yD));

	// value out of range, return upper bound
	EXPECT_EQ(7, interpolate(2.001, xD, yD));
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

TEST(PIDdigit, consistencyWithReferenceImplementation){
	// Tests the performance improved version against the default one
	unsigned long testPID = rand() % 1000000000 + 1000000000;
	for(size_t i=1; i < 8; i++)
	{
		HepPID::location loc = (HepPID::location) i;
		unsigned short newResult = HepPID::digit(loc, testPID);
		//original implementation
		int numerator = (int) std::pow(10.0,(loc-1));
		EXPECT_EQ(newResult, (HepPID::abspid(testPID)/numerator)%10);
	}
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
	Vector3d some_grid_point = origin + Vector3d(2, 3, 4) * spacing + Vector3d(spacing / 2.);
	EXPECT_EQ(some_grid_point, grid.positionFromIndex(some_index));

	grid.get(2, 3, 4) = 7;
	EXPECT_FLOAT_EQ(7., grid.getGrid()[some_index]);
	EXPECT_FLOAT_EQ(7., grid.interpolate(some_grid_point));
}

TEST(ScalarGrid, ClosestValue) {
	// Check some closest values
	ScalarGrid grid(Vector3d(0.), 2, 2, 2, 1.);
	grid.get(0, 0, 0) = 1;
	grid.get(0, 0, 1) = 2;
	grid.get(0, 1, 0) = 3;
	grid.get(0, 1, 1) = 4;
	grid.get(1, 0, 0) = 5;
	grid.get(1, 0, 1) = 6;
	grid.get(1, 1, 0) = 7;
	grid.get(1, 1, 1) = 8;

	// grid points are at 0.5 and 1.5
	EXPECT_FLOAT_EQ(1, grid.closestValue(Vector3d(0.1,  0.2, 0.6)));
	EXPECT_FLOAT_EQ(2, grid.closestValue(Vector3d(0.2, 0.1, 1.3)));
	EXPECT_FLOAT_EQ(3, grid.closestValue(Vector3d(0.3, 1.2, 0.2)));
	EXPECT_FLOAT_EQ(7, grid.closestValue(Vector3d(1.7, 1.8, 0.4)));
}

TEST(VectorGrid, Interpolation) {
	// Explicitly test trilinear interpolation
	double spacing = 2.793;
	int n = 3;
	VectorGrid grid(Vector3d(0.), n, n, n, spacing);
	grid.get(0, 0, 1) = Vector3f(1.7, 0., 0.); // set one value

	Vector3d b;

	// grid points are at (0.5, 1.5, 2.5) * spacing
	b = grid.interpolate(Vector3d(0.5, 0.5, 1.5) * spacing);
	EXPECT_FLOAT_EQ(1.7, b.x);

	b = grid.interpolate(Vector3d(0.5, 0.5, 1.4) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9, b.x);

	b = grid.interpolate(Vector3d(0.5, 0.5, 1.6) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9, b.x);

	b = grid.interpolate(Vector3d(0.5, 0.35, 1.6) * spacing);
	EXPECT_FLOAT_EQ(1.7 * 0.9 * 0.85, b.x);

	b = grid.interpolate(Vector3d(0.5, 2.65, 1.6) * spacing);
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

TEST(CylindricalProjectionMap, functions) {
	Vector3d v;
	v.setRThetaPhi(1.0, 1.2, 2.4);
	EXPECT_NEAR(v.getPhi(), 2.4, .00001);
	EXPECT_NEAR(v.getTheta(), 1.2, .000001);



	CylindricalProjectionMap cpm(24, 12);
	size_t bin = 50;
	Vector3d d = cpm.directionFromBin(bin);
	size_t bin2 = cpm.binFromDirection(d);
	EXPECT_EQ(bin, bin2);
}

TEST(EmissionMap, functions) {

	EmissionMap em(360, 180, 100, 1 * EeV, 100 * EeV);
	double e = em.energyFromBin(50);
	size_t b = em.binFromEnergy(50 * EeV);

	Vector3d d(1.0, 0.0, 0.0);

	em.fillMap(1, 50 * EeV, d);

	Vector3d d2;

	bool r = em.drawDirection(1, 50 * EeV, d2);
	EXPECT_TRUE(r);
	EXPECT_TRUE(d.getAngleTo(d2) < (2. * M_PI / 180.));

	r = em.drawDirection(1, 30 * EeV, d2);
	EXPECT_FALSE(r);

	r = em.drawDirection(2, 50 * EeV, d2);
	EXPECT_FALSE(r);

}

TEST(EmissionMap, merge) {
	EmissionMap em1, em2;
	em1.fillMap(1, 50 * EeV, Vector3d(1.0, 0.0, 0.0));
	em2.fillMap(1, 50 * EeV, Vector3d(0.0, 1.0, 0.0));
	em2.fillMap(2, 50 * EeV, Vector3d(0.0, 1.0, 0.0));

	em1.merge(&em2);

	EXPECT_EQ(em1.getMaps().size(), 2);

	ref_ptr<CylindricalProjectionMap> cpm = em1.getMap(1, 50 * EeV);
	size_t bin = cpm->binFromDirection(Vector3d(0.0, 1.0, 0.0));
	EXPECT_TRUE(cpm->getPdf()[bin] > 0);
}


TEST(Variant, copyToBuffer)
{
	double a = 23.42;
	Variant v(a);
	double b;
	v.copyToBuffer(&b);
	EXPECT_EQ(a, b);
}

TEST(Variant, stringConversion)
{
	Variant v, w;
	{
		int32_t a = 12;
		v = Variant::fromInt32(a);
		EXPECT_EQ(a, v.asInt32());

		w = Variant::fromString(v.toString(), v.getType());
		EXPECT_EQ(a, w.asInt32());
	}

	{
		int64_t a = 12;
		v = Variant::fromInt64(a);
		EXPECT_EQ(a, v.asInt64());

		w = Variant::fromString(v.toString(), v.getType());
		EXPECT_EQ(a, w.asInt64());
	}
}


TEST(Geometry, Plane)
{
	Plane p(Vector3d(0,0,1), Vector3d(0,0,1));
	EXPECT_DOUBLE_EQ(-1., p.distance(Vector3d(0, 0, 0)));
	EXPECT_DOUBLE_EQ(9., p.distance(Vector3d(1, 1, 10)));
}

TEST(Geometry, Sphere)
{
	Sphere s(Vector3d(0,0,0), 1.);
	EXPECT_DOUBLE_EQ(-1., s.distance(Vector3d(0, 0, 0)));
	EXPECT_DOUBLE_EQ(9., s.distance(Vector3d(10, 0, 0)));
}

TEST(Geometry, ParaxialBox)
{
	ParaxialBox b(Vector3d(0,0,0), Vector3d(3,4,5));
	EXPECT_NEAR(-.1, b.distance(Vector3d(0.1, 0.1, 0.1)), 1E-10);
	EXPECT_NEAR(-.1, b.distance(Vector3d(0.1, 3.8, 0.1)), 1E-10);
	EXPECT_NEAR(-.2, b.distance(Vector3d(0.9, 3.8, 0.9)), 1E-10);
	EXPECT_NEAR(7., b.distance(Vector3d(10., 0., 0.)), 1E-10);
	EXPECT_NEAR(7., b.distance(Vector3d(10., 2., 0.)), 1E-10);
	EXPECT_NEAR(8., b.distance(Vector3d(-8., 0., 0.)), 1E-10);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

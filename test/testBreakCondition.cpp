/** Unit tests for break condition, observer, boundary and tool modules */

#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/Observer.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/Tools.h"
#include "crpropa/module/RestrictToRegion.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Geometry.h"

#include "gtest/gtest.h"

namespace crpropa {

//** ========================= Break conditions ============================= */
TEST(MinimumEnergy, test) {
	MinimumEnergy minEnergy(5);
	Candidate c;

	c.current.setEnergy(5.1);
	minEnergy.process(&c);
	EXPECT_TRUE(c.isActive());

	c.current.setEnergy(4.9);
	minEnergy.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(MinimumChargeNumber, test) {
	MinimumChargeNumber minChargeNumber(20);
	Candidate c;

	c.current.setId(nucleusId(56, 26));
	minChargeNumber.process(&c);
	EXPECT_TRUE(c.isActive());
	
	c.current.setId(-nucleusId(56, 26));
	minChargeNumber.process(&c);
	EXPECT_TRUE(c.isActive());

	c.current.setId(nucleusId(4, 2));
	minChargeNumber.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
	
	c.setActive(true);
	c.removeProperty("Rejected");

	c.current.setId(-nucleusId(4, 2));
	minChargeNumber.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(MinimumEnergyPerParticleId, test) {
	MinimumEnergyPerParticleId minEnergy(1);
	minEnergy.add(22, 10);
	minEnergy.add(12, 2);
	minEnergy.add(11, 20);

	Candidate c;

	c.current.setEnergy(20);
	c.current.setId(22);
	minEnergy.process(&c);
	EXPECT_TRUE(c.isActive());

	c.current.setEnergy(5);
	c.current.setId(22);
	minEnergy.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));

	c.setActive(true);
	c.removeProperty("Rejected");

	c.current.setEnergy(10);
	c.current.setId(11);
	minEnergy.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));

	c.setActive(true);
	c.removeProperty("Rejected");

	c.current.setEnergy(5);
	c.current.setId(12);
	minEnergy.process(&c);
	EXPECT_TRUE(c.isActive());	

	c.current.setEnergy(0.1);
	c.current.setId(12);
	minEnergy.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(MaximumTrajectoryLength, test) {
	MaximumTrajectoryLength maxLength(10);
	Candidate c;

	c.setTrajectoryLength(9.9);
	maxLength.process(&c);
	EXPECT_TRUE(c.isActive());

	c.setTrajectoryLength(10.1);
	maxLength.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(MaximumTrajectoryLength, observer) {
	MaximumTrajectoryLength maxLength(12);
	maxLength.addObserverPosition(Vector3d(10, 0, 0));
	Candidate c;
	c.current.setPosition(Vector3d(5, 0, 0));

	c.setTrajectoryLength(5);
	maxLength.process(&c);
	EXPECT_TRUE(c.isActive());

	c.setTrajectoryLength(8);
	maxLength.process(&c);
	EXPECT_FALSE(c.isActive());
}

TEST(MinimumRedshift, test) {
	MinimumRedshift minZ; // default minimum redshift of 0
	Candidate c;

	c.setRedshift(0.1);
	minZ.process(&c);
	EXPECT_TRUE(c.isActive());

	c.setRedshift(0);
	minZ.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(DetectionLength, test) {
        DetectionLength detL(10);
	detL.setMakeRejectedInactive(false);
        Candidate c;
        c.current.setPosition(Vector3d(5,0,0));

        c.setTrajectoryLength(2);
        detL.process(&c);
        EXPECT_TRUE(c.isActive());
	
        c.setCurrentStep(10);
	c.setTrajectoryLength(12);
        detL.process(&c);
        EXPECT_TRUE(c.isActive());
        EXPECT_TRUE(c.hasProperty("Rejected"));
}

//** ============================= Observers ================================ */
TEST(ObserverFeature, SmallSphere) {
	// detect if the current position is inside and the previous outside of the sphere
	Observer obs;
	obs.add(new ObserverSurface(new Sphere (Vector3d(0, 0, 0), 1)));
	Candidate c;
	c.setNextStep(10);

	// no detection: particle was inside already
	c.current.setPosition(Vector3d(0.9, 0, 0));
	c.previous.setPosition(Vector3d(0.95, 0, 0));
	obs.process(&c);
	EXPECT_TRUE(c.isActive());

	// limit step
	EXPECT_NEAR(c.getNextStep(), 0.1, 0.001);

	// detection: particle just entered
	c.current.setPosition(Vector3d(0.9, 0, 0));
	c.previous.setPosition(Vector3d(1.1, 0, 0));
	obs.process(&c);
	EXPECT_FALSE(c.isActive());
}

TEST(ObserverFeature, LargeSphere) {
	// detect if the current position is outside and the previous inside of the sphere
	Observer obs;
	obs.add(new ObserverSurface(new Sphere (Vector3d(0, 0, 0), 10)));
	Candidate c;
	c.setNextStep(10);

	// no detection: particle was outside already
	c.current.setPosition(Vector3d(11, 0, 0));
	c.previous.setPosition(Vector3d(10.5, 0, 0));
	obs.process(&c);
	EXPECT_TRUE(c.isActive());

	// limit step
	EXPECT_DOUBLE_EQ(c.getNextStep(), 1);

	// detection: particle just left
	c.current.setPosition(Vector3d(11, 0, 0));
	c.previous.setPosition(Vector3d(9.5, 0, 0));
	obs.process(&c);
	EXPECT_FALSE(c.isActive());
}

TEST(ObserverFeature, Point) {
	Observer obs;
	obs.add(new Observer1D());
	Candidate c;
	c.setNextStep(10);

	// no detection, limit step
	c.current.setPosition(Vector3d(5, 0, 0));
	obs.process(&c);
	EXPECT_TRUE(c.isActive());

	// limit step
	EXPECT_DOUBLE_EQ(5, c.getNextStep());

	// detection
	c.current.setPosition(Vector3d(0, 0, 0));
	obs.process(&c);
	EXPECT_FALSE(c.isActive());
}

TEST(ObserverFeature, DetectAll) {
	// DetectAll should detect all candidates
	Observer obs;
	obs.add(new ObserverDetectAll());
	Candidate c;
	obs.process(&c);
	EXPECT_FALSE(c.isActive());
}

TEST(ObserverFeature, TimeEvolution) {
  Observer obs;
  obs.setDeactivateOnDetection(false);
  obs.setFlag("Detected", "Detected");
  //min = 5, max = min + (numb-1)*dist = 5 + 1*5 = 10, detection can happen at [5, 10]
  obs.add(new ObserverTimeEvolution(5, 5, 2));
  Candidate c;
  c.setNextStep(10);
  c.setTrajectoryLength(3);
  
  // Simulate simple detections to guarantee ObserverTimeEvolution.checkDetection is working:
  // no detection, limit next step
  obs.process(&c);
  EXPECT_TRUE(c.isActive());

  // limit step
  EXPECT_DOUBLE_EQ(2, c.getNextStep());
  
  // detection one
  c.setCurrentStep(0.1);
  c.setTrajectoryLength(5);
  obs.process(&c);
  EXPECT_TRUE(c.isActive());
  EXPECT_TRUE(c.hasProperty("Detected"));

  // no detection expected
  obs.setDeactivateOnDetection(true); // set this to true, so it deactivates if a detection happens (not expected)
  c.setTrajectoryLength(8);
  obs.process(&c);
  EXPECT_TRUE(c.isActive());

  // detection two
  c.setCurrentStep(0.1);
  c.setTrajectoryLength(10.05);
  obs.process(&c);
  EXPECT_FALSE(c.isActive());
  EXPECT_TRUE(c.hasProperty("Detected"));
}

TEST(ObserverFeature, TimeEvolutionLog) {
  Observer obs;
  obs.setDeactivateOnDetection(false);
  obs.setFlag("Detected", "Detected");
  // usage of a log scaling for the observer
  bool log = true;
  obs.add(new ObserverTimeEvolution(10, 1000, 3, log));
  Candidate c;
  // choose a stepsize that is larger then distance to next detection at 10 to check step limitation
  c.setNextStep(10);
  // set length before next detection
  c.setTrajectoryLength(3);

  // Simulate simple detections to guarantee ObserverTimeEvolution.checkDetection is working:
  // no detection, limit next step
  obs.process(&c);
  EXPECT_TRUE(c.isActive());

  // limit step (should be 10-3=7)
  EXPECT_DOUBLE_EQ(7, c.getNextStep());

  // detection one
  c.setCurrentStep(0.1);  // set small to be barely over first detection length
  c.setTrajectoryLength(10);  // set to first detection length
  obs.process(&c);
  EXPECT_TRUE(c.isActive());
  EXPECT_TRUE(c.hasProperty("Detected"));

  // no detection expected
  obs.setDeactivateOnDetection(true); // set this to true, so it deactivates if a detection happens (not expected)
  c.setTrajectoryLength(80);  // set to something between 10 and 100 (first and second detection)
  obs.process(&c);
  EXPECT_TRUE(c.isActive());
  obs.setDeactivateOnDetection(false); // reset to false again for future detection

  // detection two
  c.setCurrentStep(0.1);
  c.setTrajectoryLength(100);
  obs.process(&c);
  EXPECT_TRUE(c.isActive());
  EXPECT_TRUE(c.hasProperty("Detected"));

  // detection two
  obs.setDeactivateOnDetection(true);  // deactivate here since it is the last detection
  c.setCurrentStep(0.1);
  c.setTrajectoryLength(1000.05);
  obs.process(&c);
  EXPECT_FALSE(c.isActive());  // not active anymore
  EXPECT_TRUE(c.hasProperty("Detected"));
}

TEST(ObserverFeature, TimeEvolutionArray) {
  // here it should be tested if the observer can be constructed with an array
  std::vector<double> times = {1, 2, 3}; 
  ObserverTimeEvolution obs(times);
  EXPECT_FALSE(obs.empty());
  EXPECT_TRUE(times == obs.getTimes());  // element wise comparison

  times.push_back(4);
  obs.addTime(4);
  EXPECT_FALSE(obs.empty());
  EXPECT_TRUE(times == obs.getTimes());

  // test clear:
  obs.clear();
  EXPECT_TRUE(obs.empty());

  // test addTimeRange for linear ranges
  times = {5, 6, 7, 8, 9, 10};
  obs.clear();  // empty detList
  EXPECT_TRUE(obs.empty());
  obs.addTimeRange(5, 10, 6, false);
  EXPECT_FALSE(obs.empty());
  EXPECT_TRUE(times == obs.getTimes());

  // test addTimeRange for logarithmic ranges
  times = {1, 10, 100, 1000};
  obs.clear();  // empty detList
  EXPECT_TRUE(obs.empty());
  obs.addTimeRange(1, 1000, 4, true);
  EXPECT_FALSE(obs.empty());
  // should be equal to above times array, but isnt, even though the values are the same
  for (int i=0; i<times.size(); i++)
    EXPECT_NEAR(times[i], obs.getTimes()[i], 0.01);

  // now check if constructDetListIfEmpty is working properly:
  ObserverTimeEvolution obs2(5, 10, 6, false);
  times = {5, 6, 7, 8, 9, 10};
  
  // check if no array is created while calling getTimes
  EXPECT_TRUE(obs2.empty());
  EXPECT_TRUE(times == obs2.getTimes());
  EXPECT_TRUE(obs2.empty());

  // check if array is created when calling addTime without array
  times.push_back(11);
  obs2.addTime(11);
  EXPECT_FALSE(obs2.empty());
  EXPECT_TRUE(times == obs2.getTimes());
}

//** ========================= Boundaries =================================== */
TEST(PeriodicBox, high) {
	// Tests if the periodical boundaries place the particle back inside the box and translate the initial position accordingly.
	Vector3d origin(2, 2, 2);
	Vector3d size(2, 2, 2);
	PeriodicBox box(origin, size);

	Candidate c;
	c.current.setPosition(Vector3d(4.5, 4.3, 4.4));
	c.created.setPosition(Vector3d(3, 3, 3));

	box.process(&c);

	EXPECT_DOUBLE_EQ(2.5, c.current.getPosition().x);
	EXPECT_DOUBLE_EQ(1, c.created.getPosition().x);
	EXPECT_DOUBLE_EQ(2.3, c.current.getPosition().y);
	EXPECT_DOUBLE_EQ(1, c.created.getPosition().y);
	EXPECT_DOUBLE_EQ(2.4, c.current.getPosition().z);
	EXPECT_DOUBLE_EQ(1, c.created.getPosition().z);
}

TEST(PeriodicBox, low) {
	// Tests if the periodical boundaries place the particle back inside the box and translate the initial position accordingly.
	Vector3d origin(0, 0, 0);
	Vector3d size(2, 2, 2);
	PeriodicBox box(origin, size);

	Candidate c;
	c.current.setPosition(Vector3d(-2.5, -0.3, -0.4));
	c.created.setPosition(Vector3d(1, 1, 1));

	box.process(&c);

	EXPECT_DOUBLE_EQ(1.5, c.current.getPosition().x);
	EXPECT_DOUBLE_EQ(5, c.created.getPosition().x);
	EXPECT_DOUBLE_EQ(1.7, c.current.getPosition().y);
	EXPECT_DOUBLE_EQ(3, c.created.getPosition().y);
	EXPECT_DOUBLE_EQ(1.6, c.current.getPosition().z);
	EXPECT_DOUBLE_EQ(3, c.created.getPosition().z);
}

TEST(ReflectiveShell, inside) {
	// Tests if the reflective boundaries place the particle back inside the shell
	Vector3d center(0, 0, 0);
	double radius = 100;
	ReflectiveShell shell(center, radius);

	Candidate c;
	c.setCurrentStep(20);
	c.previous.setPosition(Vector3d(80, 20, 30));
	c.previous.setDirection(Vector3d(10, -1, -1));
	// un-reflected new position (outside shell) after full step
	Vector3d currentPosition = c.previous.getPosition() + c.previous.getDirection() * c.getCurrentStep() / c.previous.getDirection().getR();
	c.current.setPosition(currentPosition);
	c.current.setDirection(Vector3d(10, -1, -1));
	// process reflection
	shell.process(&c);

	// expected position & direction after reflection
	// calculated by hand with the same algorithm as implemented in Boundary.cpp
	EXPECT_NEAR(90.06356247, c.current.getPosition().x, 1e-7);
	EXPECT_NEAR(16.09256317, c.current.getPosition().y, 1e-7);
	EXPECT_NEAR(25.05646296, c.current.getPosition().z, 1e-7);

	EXPECT_NEAR(-0.67179589, c.current.getDirection().x, 1e-7);
	EXPECT_NEAR(-0.42786503, c.current.getDirection().y, 1e-7);
	EXPECT_NEAR(-0.60466668, c.current.getDirection().z, 1e-7);
}

TEST(ReflectiveBox, high) {
	// Tests if the reflective boundaries place the particle back inside the box and translate the initial position accordingly.
	// Also the initial and final directions are to be reflected
	Vector3d origin(10, 10, 10);
	Vector3d size(10, 20, 20);
	ReflectiveBox box(origin, size);

	Candidate c;
	c.source.setPosition(Vector3d(16, 17, 18));
	c.source.setDirection(Vector3d(1, 1.6, 1.8));
	c.created.setPosition(Vector3d(15, 15, 15));
	c.created.setDirection(Vector3d(0, 0.6, 0.8));
	c.previous.setPosition(Vector3d(15, 15, 29.5));
	c.previous.setDirection(Vector3d(0, 0.6, 0.8));
	c.current.setPosition(Vector3d(15, 15, 30.5));
	c.current.setDirection(Vector3d(0, 0.6, 0.8));

	box.process(&c);

	EXPECT_DOUBLE_EQ(16, c.source.getPosition().x);
	EXPECT_DOUBLE_EQ(17, c.source.getPosition().y);
	EXPECT_DOUBLE_EQ(42, c.source.getPosition().z);

	EXPECT_DOUBLE_EQ(15, c.created.getPosition().x);
	EXPECT_DOUBLE_EQ(15, c.created.getPosition().y);
	EXPECT_DOUBLE_EQ(45, c.created.getPosition().z);

	EXPECT_DOUBLE_EQ(15, c.previous.getPosition().x);
	EXPECT_DOUBLE_EQ(15, c.previous.getPosition().y);
	EXPECT_DOUBLE_EQ(30.5, c.previous.getPosition().z);

	EXPECT_DOUBLE_EQ(15, c.current.getPosition().x);
	EXPECT_DOUBLE_EQ(15, c.current.getPosition().y);
	EXPECT_DOUBLE_EQ(29.5, c.current.getPosition().z);

	EXPECT_DOUBLE_EQ(0, c.created.getDirection().x);
	EXPECT_DOUBLE_EQ(0.6, c.created.getDirection().y);
	EXPECT_DOUBLE_EQ(-0.8, c.created.getDirection().z);

	EXPECT_DOUBLE_EQ(0, c.previous.getDirection().x);
	EXPECT_DOUBLE_EQ(0.6, c.previous.getDirection().y);
	EXPECT_DOUBLE_EQ(-0.8, c.previous.getDirection().z);

	EXPECT_DOUBLE_EQ(0, c.current.getDirection().x);
	EXPECT_DOUBLE_EQ(0.6, c.current.getDirection().y);
	EXPECT_DOUBLE_EQ(-0.8, c.current.getDirection().z);
}

TEST(CubicBoundary, inside) {
	CubicBoundary cube(Vector3d(0, 0, 0), 10);
	Candidate c;
	c.current.setPosition(Vector3d(9, 5, 5));
	cube.process(&c);
	EXPECT_TRUE(c.isActive());
}

TEST(CubicBoundary, outside) {
	CubicBoundary cube(Vector3d(0, 0, 0), 10);
	Candidate c;
	c.current.setPosition(Vector3d(10.1, 5, 5));
	cube.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(CubicBoundary, limitStepLower) {
	CubicBoundary cube(Vector3d(10, 10, 10), 10);
	cube.setLimitStep(true);
	cube.setMargin(1);
	Candidate c;
	c.current.setPosition(Vector3d(15, 15, 10.5));
	c.setNextStep(100);
	cube.process(&c);
	EXPECT_DOUBLE_EQ(1.5, c.getNextStep());
}

TEST(CubicBoundary, limitStepUpper) {
	CubicBoundary cube(Vector3d(-10, -10, -10), 10);
	cube.setLimitStep(true);
	cube.setMargin(1);
	Candidate c;
	c.current.setPosition(Vector3d(-5, -5, -0.5));
	c.setNextStep(100);
	cube.process(&c);
	EXPECT_DOUBLE_EQ(1.5, c.getNextStep());
}

TEST(SphericalBoundary, inside) {
	SphericalBoundary sphere(Vector3d(0, 0, 0), 10);
	Candidate c;
	c.current.setPosition(Vector3d(9, 0, 0));
	sphere.process(&c);
	EXPECT_TRUE(c.isActive());
	EXPECT_FALSE(c.hasProperty("Rejected"));
}

TEST(SphericalBoundary, outside) {
	SphericalBoundary sphere(Vector3d(0, 0, 0), 10);
	sphere.setRejectFlag("I passed the galactic border", "Nothing happened");
	Candidate c;
	c.current.setPosition(Vector3d(0, -10.1, 0));
	sphere.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("I passed the galactic border"));
}

TEST(SphericalBoundary, limitStep) {
	SphericalBoundary sphere(Vector3d(0, 0, 0), 10);
	sphere.setLimitStep(true);
	sphere.setMargin(1);
	Candidate c;
	c.setNextStep(100);
	c.current.setPosition(Vector3d(0, 0, 9.5));
	sphere.process(&c);
	EXPECT_DOUBLE_EQ(1.5, c.getNextStep());
}

TEST(EllipsoidalBoundary, inside) {
	EllipsoidalBoundary ellipsoid(Vector3d(-5, 0, 0), Vector3d(5, 0, 0), 15);
	Candidate c;
	c.current.setPosition(Vector3d(3, 2, 0));
	ellipsoid.process(&c);
	EXPECT_TRUE(c.isActive());
	EXPECT_FALSE(c.hasProperty("Rejected"));
}

TEST(EllipsoidalBoundary, outside) {
	EllipsoidalBoundary ellipsoid(Vector3d(-5, 0, 0), Vector3d(5, 0, 0), 15);
	Candidate c;
	c.current.setPosition(Vector3d(0, 25, 0));
	ellipsoid.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(EllipsoidalBoundary, limitStep) {
	EllipsoidalBoundary ellipsoid(Vector3d(-5, 0, 0), Vector3d(5, 0, 0), 15);
	ellipsoid.setLimitStep(true);
	ellipsoid.setMargin(0.5);
	Candidate c;
	c.setNextStep(2);
	c.current.setPosition(Vector3d(7, 0, 0));
	ellipsoid.process(&c);
	EXPECT_DOUBLE_EQ(c.getNextStep(), 1.5);
}

TEST(CylindricalBoundary, inside) {
        CylindricalBoundary cylinder(Vector3d(0, 0, 0), 2, 15);
	Candidate c;
	c.current.setPosition(Vector3d(6, -3, 0.5));
	cylinder.process(&c);
	EXPECT_TRUE(c.isActive());
	EXPECT_FALSE(c.hasProperty("Rejected"));
}

TEST(CylindricalBoundary, outside) {
        CylindricalBoundary cylinder(Vector3d(0, 0, 0), 2, 15);
	Candidate c;
	c.current.setPosition(Vector3d(6, -3, 1.5));
	cylinder.process(&c);
	EXPECT_FALSE(c.isActive());
	EXPECT_TRUE(c.hasProperty("Rejected"));
}

TEST(CylindricalBoundary, limitStep) {
        CylindricalBoundary cylinder(Vector3d(0, 0, 0), 2, 15);
	cylinder.setLimitStep(true);
	cylinder.setMargin(0.5);
	Candidate c;
	c.setNextStep(2);
	c.current.setPosition(Vector3d(7, 0, 0));
	cylinder.process(&c);
	EXPECT_DOUBLE_EQ(c.getNextStep(), 1.5);
}

TEST(RestrictToRegion, RestrictToRegion) {

	ref_ptr<Observer> obs = new Observer();
	obs->add(new ObserverDetectAll());
	RestrictToRegion R(obs, new Sphere(Vector3d(0, 0, 0), 10));

	Candidate c;
	c.previous.setPosition(Vector3d(13,0,0));
	c.current.setPosition(Vector3d(12,0,0));
	R.process(&c);
	EXPECT_TRUE(c.isActive());
	c.current.setPosition(Vector3d(9,0,0));
	R.process(&c);
	EXPECT_FALSE(c.isActive());
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

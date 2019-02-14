#include "crpropa/Candidate.h"
#include "crpropa/ParticleID.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/PropagationBP.h"
#include "crpropa/module/PropagationCK.h"

#include "gtest/gtest.h"

#include <string>
#include <iostream>




namespace crpropa {

TEST(testSimplePropagation, step) {
	double minStep = 20;
	double maxStep = 100;
	SimplePropagation propa(minStep, maxStep);

	ParticleState p;
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	c.setNextStep(10);

	propa.process(&c);

	EXPECT_EQ(minStep, c.getCurrentStep());
	EXPECT_EQ(maxStep, c.getNextStep());
	EXPECT_EQ(Vector3d(0,  0, 0), c.created.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.created.getDirection());
	EXPECT_EQ(Vector3d(0,  0, 0), c.previous.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.previous.getDirection());
	EXPECT_EQ(Vector3d(0, 20, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.current.getDirection());
}


TEST(testPropagationCK, zeroField) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 0)));

	double minStep = 0.1 * kpc;
	propa.setMinimumStep(minStep);

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * minStep, c.getNextStep());  // acceleration by factor 5
}


TEST(testPropagationCK, proton) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	double minStep = 0.1 * kpc;
	propa.setMinimumStep(minStep);

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * minStep, c.getNextStep());  // acceleration by factor 5
}

TEST(testPropagationCK, neutron) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
	propa.setMinimumStep(1 * kpc);
	propa.setMaximumStep(42 * Mpc);

	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(1 * kpc, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(42 * Mpc, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}

TEST(testPropagationBP, zeroField) {
	PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 0)), 1 * kpc);

	double minStep = 0.1 * kpc;
	propa.setMinimumStep(minStep);
	propa.setTolerance(0.42);

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * minStep, c.getNextStep());  // acceleration by factor 5

}

TEST(testPropagationBP, exceptions)
{

    // minStep should be smaller than maxStep
    EXPECT_THROW(PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 10 , 0), std::runtime_error);
    // Too large tolerance: tolerance should be between 0 and 1
    EXPECT_THROW(PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 10 * kpc , 20 * kpc, 42.), std::runtime_error);

    PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

    // set maximum step, so that it can be tested what happens if a larger minStep is set.
    propa.setMaximumStep(1 * Mpc);

    // this tests _that_ the expected exception is thrown
    EXPECT_THROW(propa.setTolerance(2.), std::runtime_error);
    EXPECT_THROW(propa.setMinimumStep(-1.), std::runtime_error);
    EXPECT_THROW(propa.setMinimumStep(2 * Mpc), std::runtime_error);

    // set minimum step, so that it can be tested what happens if a smaller maxStep is set.
    propa.setMinimumStep(0.5 * Mpc);

    EXPECT_THROW(propa.setMaximumStep(0.1 * Mpc), std::runtime_error);
}


TEST(testPropagationBP, constructor) {
    // Test construction and parameters
    ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(0, 0, 1 * nG));

    double minStep = 1.;
    double maxStep = 100.;
    double tolerance = 0.01;

    PropagationBP propa(bField, minStep, maxStep, tolerance);

    EXPECT_EQ(minStep, propa.getMinimumStep());
    EXPECT_EQ(maxStep, propa.getMaximumStep());
    EXPECT_EQ(tolerance, propa.getTolerance());
    EXPECT_EQ(bField, propa.getField());

    // Update parameters
    minStep = 10.;
    maxStep = 10.;
    propa.setTolerance(0.0001);
    bField = new UniformMagneticField(Vector3d(10 * nG, 0, 1 * nG));

    propa.setTolerance(tolerance);
    propa.setMinimumStep(minStep);
    propa.setMaximumStep(maxStep);
    propa.setField(bField);

    EXPECT_EQ(minStep, propa.getMinimumStep());
    EXPECT_EQ(maxStep, propa.getMaximumStep());
    EXPECT_EQ(tolerance, propa.getTolerance());
    EXPECT_EQ(bField, propa.getField());

    // Test the fixed step size version of the Boris push
    minStep = 10. * kpc;
    PropagationBP propaBP(bField, minStep);
    EXPECT_EQ(propaBP.getMaximumStep(), propaBP.getMaximumStep());

    // The propagation should be initialized with the default constructor
    PropagationBP propaBPField(bField);
    EXPECT_EQ(propaBPField.getMaximumStep(), propaBPField.getMaximumStep());
    EXPECT_EQ(propaBPField.getMaximumStep(), 1 * kpc);
}

// Test if the step size is reduced correctly if the error is too large with respect to the tolerance: r > 1
TEST(testPropagationBP, reduceStep) {
    PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 100 * nG)), 1 * kpc);

    double minStep = 0.1 * kpc;
    double maxStep = 1 * Gpc;
    propa.setMinimumStep(minStep);
    propa.setMaximumStep(maxStep);
    // small tolerance leads to large values of r
    propa.setTolerance(1e-15);

    ParticleState p;
    p.setId(nucleusId(1, 1));
    p.setEnergy(100 * TeV);
    p.setPosition(Vector3d(0, 0, 0));
    p.setDirection(Vector3d(0, 1, 0));
    Candidate c(p);
    // large step leads to large errors and thus in combination with the low tolerance to high values of r
    c.setNextStep(maxStep);

    propa.process(&c);

    // adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
    EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
    EXPECT_DOUBLE_EQ(minStep, c.getNextStep());  // stay at minimum step because of large r due to small tolerance

}

// Test if the step size is increased correctly if the error is small with respect to the tolerance: r < 1
TEST(testPropagationBP, increaseStep) {
    PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 1 * kpc);

    double minStep = 0.001 * pc;
    double maxStep = 3.125 * pc;
    propa.setMinimumStep(minStep);
    propa.setMaximumStep(maxStep);
    // large tolerance leads to small values of r. Consequently, the step size can be increased.
    propa.setTolerance(0.9);

    ParticleState p;
    p.setId(nucleusId(1, 1));
    p.setEnergy(100 * EeV);
    p.setPosition(Vector3d(0, 0, 0));
    p.setDirection(Vector3d(0, 1, 0));
    Candidate c(p);

    // each step the step size can be increased by a factor of 5.
    for (int i = 1; i < 6; i++){
        propa.process(&c);
        EXPECT_DOUBLE_EQ(minStep*pow(5, i) / pc, c.getNextStep()/pc);
    }
    // after 5 steps the maxStep is reached. The current step is, however, less.
    EXPECT_DOUBLE_EQ(maxStep/pc/5., c.getCurrentStep()/pc);
    EXPECT_DOUBLE_EQ(maxStep/pc, c.getNextStep()/pc);
}

TEST(testPropagationBP, proton) {
	PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	double step = 0.01 * kpc;
	propa.setMinimumStep(step);
	propa.setMaximumStep(10*step);
	propa.setTolerance(0.00001);

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(step, c.getCurrentStep());  // perform step
	EXPECT_DOUBLE_EQ(5 * step, c.getNextStep());  // acceleration by factor 5
}


// Test the numerical results for parallel magnetic field lines along the z-axis
TEST(testPropagationBP, gyration) {
    PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

    double step = 10. * Mpc;  // gyroradius is 108.1 Mpc
    propa.setMaximumStep(step);
    propa.setMinimumStep(step);


    ParticleState p;
    p.setId(nucleusId(1, 1));
    p.setEnergy(100 * EeV);
    p.setPosition(Vector3d(0, 0, 0));
    p.setDirection(Vector3d(1, 1, 1));
    Candidate c(p);
    c.setNextStep(0);
    propa.process(&c);

    double dirX = c.current.getDirection().x;
    double dirY = c.current.getDirection().y;
    double dirZ = c.current.getDirection().z;
    double posZ = c.current.getPosition().z;

    // Test if the analytical solution is achieved of the components of the momentum with the Boris push as expected in
    // the background magnetic field.
    EXPECT_DOUBLE_EQ(2 / 3., dirX * dirX + dirY * dirY);  // constant momentum in the perpendicular plane to background magnetic field field
    EXPECT_DOUBLE_EQ(1 / 3., dirZ * dirZ);  // constant momentum parallel to the background magnetic field
    EXPECT_DOUBLE_EQ( step * step / 3., posZ * posZ);  // constant velocity parallel to the background magnetic field

    // Nine new steps to have finally propagated the particle ten times
    for (int i = 0; i < 9; i++){
        propa.process(&c);
    }

    dirX = c.current.getDirection().x;
    dirY = c.current.getDirection().y;
    dirZ = c.current.getDirection().z;
    posZ = c.current.getPosition().z;

    // Compare the numerical solutions after ten steps with the analytical solution of the trajectories
    EXPECT_DOUBLE_EQ(2 / 3., dirX * dirX + dirY * dirY);  // constant momentum in the perpendicular plane to background magnetic field field
    EXPECT_DOUBLE_EQ(1 / 3., dirZ * dirZ);  // constant momentum parallel to the background magnetic field
    EXPECT_DOUBLE_EQ(100 * step * step / 3., posZ * posZ);  // constant velocity parallel to the background magnetic field
}


TEST(testPropagationBP, neutron) {
	PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	propa.setMinimumStep(1 * kpc);
	propa.setMaximumStep(1 * kpc);

	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(1 * kpc, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(1 * kpc, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

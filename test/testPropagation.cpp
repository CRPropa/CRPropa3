#include "crpropa/Candidate.h"
#include "crpropa/ParticleID.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/PropagationBP.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/magneticField/turbulentField/PlaneWaveTurbulence.h"

#include "gtest/gtest.h"

#include <string>
#include <iostream>

namespace crpropa {

TEST(testSimplePropagation, timeStep) {
	ParticleState p;
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);

	double minStep = 20;  // [s]
	double maxStep = 100;  // [s]
	SimplePropagation propa(true, minStep, maxStep);
	c.setNextStep(10);

	propa.process(&c);

	EXPECT_EQ(minStep, c.getCurrentStep());
	EXPECT_EQ(maxStep, c.getNextStep());
	EXPECT_EQ(Vector3d(0,  0, 0), c.created.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.created.getDirection());
	EXPECT_EQ(Vector3d(0,  0, 0), c.previous.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.previous.getDirection());
	EXPECT_EQ(Vector3d(0, 20*p.getVelocity().getR(), 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.current.getDirection());
}


TEST(testSimplePropagation, step) {
	ParticleState p;
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);

	double minStep = 20;  // [m]
	double maxStep = 100;  // [m]
	SimplePropagation propa(minStep, maxStep);
	c.setNextStep(10/c_light);

	propa.process(&c);

	EXPECT_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());
	EXPECT_EQ(propa.getMaximumTimeStep(), c.getNextStep());
	EXPECT_EQ(Vector3d(0,  0, 0), c.created.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.created.getDirection());
	EXPECT_EQ(Vector3d(0,  0, 0), c.previous.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.previous.getDirection());
	EXPECT_EQ(Vector3d(0, 20, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0,  1, 0), c.current.getDirection());
}


TEST(testPropagationCK, zeroFieldTimeStep) {
	// test with
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 0)));
	
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);
	
	double minStep = 0.1 * kpc/c.getVelocity();
	propa.setMinimumTimeStep(minStep);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * minStep, c.getNextStep());  // acceleration by factor 5
}


TEST(testPropagationCK, zeroField) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 0)));
	
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);
	
	double minStep = 0.1 * kpc;
	propa.setMinimumStep(minStep);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * propa.getMinimumTimeStep(), c.getNextStep());  // acceleration by factor 5
}


#ifndef CRPROPA_TESTS_SKIP_EXCEPTIONS
TEST(testPropagationCK, exceptionsTimeStep) {
	// minStep should be smaller than maxStep
	EXPECT_THROW(PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 0.42, 10 , 0), std::runtime_error);
	// Too large tolerance: tolerance should be between 0 and 1
	EXPECT_THROW(PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 42., 10 * kiloyear , 20 * kiloyear), std::runtime_error);

	PropagationCK propa(1, 1 * Megayear, 1 * Megayear,  new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	// set maximum step, so that it can be tested what happens if a larger minStep is set.
	propa.setMaximumTimeStep(1 * Megayear);

	// this tests _that_ the expected exception is thrown
	EXPECT_THROW(propa.setTolerance(2.), std::runtime_error);
	EXPECT_THROW(propa.setMinimumTimeStep(-1.), std::runtime_error);
	EXPECT_THROW(propa.setMinimumTimeStep(2 * Megayear), std::runtime_error);

	// set minimum step, so that it can be tested what happens if a smaller maxStep is set.
	propa.setMinimumTimeStep(0.5 * Megayear);

	EXPECT_THROW(propa.setMaximumTimeStep(0.1 * Megayear), std::runtime_error);
}


TEST(testPropagationCK, exceptions) {
	// minStep should be smaller than maxStep
	EXPECT_THROW(PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 0.42, 10 , 0), std::runtime_error);
	// Too large tolerance: tolerance should be between 0 and 1
	EXPECT_THROW(PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 42., 10 * kpc , 20 * kpc), std::runtime_error);

	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

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
#endif


TEST(testPropagationCK, constructorTimeStep) {
	// Test construction and parameters
	// default parameters:
	ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(0, 0, 1 * nG));
	double minStep = 1.;
	double maxStep = 100.;
	double tolerance = 0.01;

	{  // check default constructor (only supports length steps)
		PropagationCK propa;

		EXPECT_FALSE(propa.getField().valid());
		EXPECT_EQ(1 * kiloparsec, propa.getMinimumStep());
		EXPECT_EQ(1 * gigaparsec, propa.getMaximumStep());
		EXPECT_EQ(1e-4, propa.getTolerance());
	}

	{  // test fixed stepsize constructor
		PropagationCK propa(0.01, minStep, minStep, bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(minStep, propa.getMaximumTimeStep());
		EXPECT_EQ(0.01, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test construction with minStep and maxStep
		PropagationCK propa(tolerance, minStep, maxStep, bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(maxStep, propa.getMaximumTimeStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test update parameters
		PropagationCK propa;

		minStep = 10.;
		maxStep = 10.;
		tolerance = 0.0001;
		bField = new UniformMagneticField(Vector3d(10 * nG, 0, 1 * nG));

		propa.setTolerance(tolerance);
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		propa.setField(bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(maxStep, propa.getMaximumTimeStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}
}


TEST(testPropagationCK, constructor) {
	// Test construction and parameters
	// default parameters:
	ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(0, 0, 1 * nG));
	double minStep = 1.;
	double maxStep = 100.;
	double tolerance = 0.01;

	{  // check default constructor (only supports length steps)
		PropagationCK propa;

		EXPECT_FALSE(propa.getField().valid());
		EXPECT_EQ(1 * kiloparsec, propa.getMinimumStep());
		EXPECT_EQ(1 * gigaparsec, propa.getMaximumStep());
		EXPECT_EQ(1e-4, propa.getTolerance());
	}

	{  // test fixed stepsize constructor
		PropagationCK propa(bField, 0.01, minStep, minStep);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(minStep, propa.getMaximumStep());
		EXPECT_EQ(0.01, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test construction with minStep and maxStep
		PropagationCK propa(bField, tolerance, minStep, maxStep);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(maxStep, propa.getMaximumStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test update parameters
		PropagationCK propa;

		minStep = 10.;
		maxStep = 10.;
		tolerance = 0.0001;
		bField = new UniformMagneticField(Vector3d(10 * nG, 0, 1 * nG));

		propa.setTolerance(tolerance);
		propa.setMinimumStep(minStep);
		propa.setMaximumStep(maxStep);
		propa.setField(bField);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(maxStep, propa.getMaximumStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}
}


// Test if the step size is reduced correctly if the error is too large with respect to the tolerance: r > 1
TEST(testPropagationCK, reduceTimeStep) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 100 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * TeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);

	double minStep = 0.1 * kpc/c.getVelocity();
	double maxStep = 1 * Gpc/c.getVelocity();
	propa.setMinimumTimeStep(minStep);
	propa.setMaximumTimeStep(maxStep);
	// small tolerance leads to large values of r
	propa.setTolerance(1e-15);
	// large step leads to large errors and thus in combination with the low tolerance to high values of r
	c.setNextStep(propa.getMaximumTimeStep());

	propa.process(&c);

	// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
}


// Test if the step size is reduced correctly if the error is too large with respect to the tolerance: r > 1
TEST(testPropagationCK, reduceStep) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 100 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * TeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);

	double minStep = 0.1 * kpc;
	double maxStep = 1 * Gpc;
	propa.setMinimumStep(minStep);
	propa.setMaximumStep(maxStep);
	// small tolerance leads to large values of r
	propa.setTolerance(1e-15);
	// large step leads to large errors and thus in combination with the low tolerance to high values of r
	c.setNextStep(maxStep);

	propa.process(&c);

	// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
}


// Test if the step size is increased correctly if the error is small with respect to the tolerance: r < 1
TEST(testPropagationCK, increaseTimeStep) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	
	double minStep = 0.001 * pc/c.getVelocity();
	double maxStep = 3.125 * pc/c.getVelocity();
	propa.setMinimumTimeStep(minStep);
	propa.setMaximumTimeStep(maxStep);
	// large tolerance leads to small values of r. Consequently, the step size can be increased.
	propa.setTolerance(0.9);

	// each step the step size can be increased by a factor of 5.
	for (int i = 1; i < 6; i++){
		propa.process(&c);
		EXPECT_DOUBLE_EQ(minStep*pow(5, i), c.getNextStep());
	}
	// after 5 steps the maxStep is reached. The current step is, however, less.
	EXPECT_DOUBLE_EQ(maxStep/5., c.getCurrentStep());
	EXPECT_DOUBLE_EQ(maxStep, c.getNextStep());
}


// Test if the step size is increased correctly if the error is small with respect to the tolerance: r < 1
TEST(testPropagationCK, increaseStep) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	
	double minStep = 0.001 * pc;
	double maxStep = 3.125 * pc;
	propa.setMinimumStep(minStep);
	propa.setMaximumStep(maxStep);
	// large tolerance leads to small values of r. Consequently, the step size can be increased.
	propa.setTolerance(0.9);

	// each step the step size can be increased by a factor of 5.
	for (int i = 1; i < 6; i++){
		propa.process(&c);
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep()*pow(5, i), c.getNextStep());
	}
	// after 5 steps the maxStep is reached. The current step is, however, less.
	EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep()/5., c.getCurrentStep());
	EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep(), c.getNextStep());
}


TEST(testPropagationCK, protonTimeStep) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);
	
	double minStep = 0.1 * kpc/c.getVelocity();
	propa.setMinimumTimeStep(minStep);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * minStep, c.getNextStep());  // acceleration by factor 5
}


TEST(testPropagationCK, proton) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);
	
	double minStep = 0.1 * kpc;
	propa.setMinimumStep(minStep);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * propa.getMinimumTimeStep(), c.getNextStep());  // acceleration by factor 5
}


// Test the numerical results for parallel magnetic field lines along the z-axis
TEST(testPropagationCK, gyrationTimeStep) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(1, 1, 1));
	Candidate c(p);
	c.setNextStep(0);
	
	double step = 10. * Mpc/c.getVelocity();  // gyroradius is 108.1 Mpc
	propa.setMaximumTimeStep(step);
	propa.setMinimumTimeStep(step);

	propa.process(&c);

	double dirX = c.current.getDirection().x;
	double dirY = c.current.getDirection().y;
	double dirZ = c.current.getDirection().z;
	double posZ = c.current.getPosition().z;

	// Test if the analytical solution is achieved for the components of the momentum with the CK method as expected in
	// the background magnetic field.
	double precision = 1e-7;
	double expected = 2 / 3.;
	EXPECT_NEAR(expected, dirX * dirX + dirY * dirY, expected * precision);  // constant momentum in the plane perpendicular to background magnetic field field
	expected = 1 / 3.;
	EXPECT_NEAR(expected, dirZ * dirZ, expected * precision);  // constant momentum parallel to the background magnetic field
	expected = step * step * c.getVelocity() * c.getVelocity() / 3.;
	EXPECT_NEAR(expected, posZ * posZ, expected * precision);  // constant velocity parallel to the background magnetic field

	// Nine new steps to have finally propagated the particle ten times
	for (int i = 0; i < 9; i++){
		propa.process(&c);
	}

	dirX = c.current.getDirection().x;
	dirY = c.current.getDirection().y;
	dirZ = c.current.getDirection().z;
	posZ = c.current.getPosition().z;

	// Compare the numerical solutions after ten steps with the analytical solution of the trajectories
	expected = 2 / 3.;
	EXPECT_NEAR(expected, dirX * dirX + dirY * dirY, expected * precision);  // constant momentum in the plane perpendicular to background magnetic field field
	expected = 1 / 3.;
	EXPECT_NEAR(expected, dirZ * dirZ, expected * precision);  // constant momentum parallel to the background magnetic field
	expected = 100 * step * step * c.getVelocity() * c.getVelocity() / 3.;
	EXPECT_NEAR(expected, posZ * posZ, expected * precision);  // constant velocity parallel to the background magnetic field
}


// Test the numerical results for parallel magnetic field lines along the z-axis
TEST(testPropagationCK, gyration) {
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(1, 1, 1));
	Candidate c(p);
	c.setNextStep(0);
	
	double step = 10. * Mpc;  // gyroradius is 108.1 Mpc
	propa.setMaximumStep(step);
	propa.setMinimumStep(step);

	propa.process(&c);

	double dirX = c.current.getDirection().x;
	double dirY = c.current.getDirection().y;
	double dirZ = c.current.getDirection().z;
	double posZ = c.current.getPosition().z;

	// Test if the analytical solution is achieved for the components of the momentum with the CK method as expected in
	// the background magnetic field.
	double precision = 1e-7;
	double expected = 2 / 3.;
	EXPECT_NEAR(expected, dirX * dirX + dirY * dirY, expected * precision);  // constant momentum in the plane perpendicular to background magnetic field field
	expected = 1 / 3.;
	EXPECT_NEAR(expected, dirZ * dirZ, expected * precision);  // constant momentum parallel to the background magnetic field
	expected = step * step / 3.;
	EXPECT_NEAR(expected, posZ * posZ, expected * precision);  // constant velocity parallel to the background magnetic field

	// Nine new steps to have finally propagated the particle ten times
	for (int i = 0; i < 9; i++){
		propa.process(&c);
	}

	dirX = c.current.getDirection().x;
	dirY = c.current.getDirection().y;
	dirZ = c.current.getDirection().z;
	posZ = c.current.getPosition().z;

	// Compare the numerical solutions after ten steps with the analytical solution of the trajectories
	expected = 2 / 3.;
	EXPECT_NEAR(expected, dirX * dirX + dirY * dirY, expected * precision);  // constant momentum in the plane perpendicular to background magnetic field field
	expected = 1 / 3.;
	EXPECT_NEAR(expected, dirZ * dirZ, expected * precision);  // constant momentum parallel to the background magnetic field
	expected = 100 * step * step / 3.;
	EXPECT_NEAR(expected, posZ * posZ, expected * precision);  // constant velocity parallel to the background magnetic field
}


TEST(testPropagationCK, neutronTimeStep) {
	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
	double minStep = 1*kpc/c.getVelocity();
	double maxStep = 42*Mpc/c.getVelocity();
	propa.setMinimumTimeStep(minStep);
	propa.setMaximumTimeStep(maxStep);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(maxStep, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}


TEST(testPropagationCK, neutron) {
	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	
	PropagationCK propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
	double minStep = 1*kpc;
	double maxStep = 42*Mpc;
	propa.setMinimumStep(minStep);
	propa.setMaximumStep(maxStep);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());
	EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep(), c.getNextStep());
	EXPECT_EQ(Vector3d(0, 1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}


TEST(testPropagationBP, zeroFieldTimeStep) {
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	PropagationBP propa(
		1. * kiloyear, 
		new UniformMagneticField(Vector3d(0, 0, 0))
	);
	
	double minStep = 0.01 * kpc/c.getVelocity();
	propa.setMinimumTimeStep(minStep);
	propa.setTolerance(0.42);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(minStep, c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * minStep, c.getNextStep());  // acceleration by factor 5
}


TEST(testPropagationBP, zeroField) {
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	Candidate c(p);
	c.setNextStep(0);

	PropagationBP propa(
		new UniformMagneticField(Vector3d(0, 0, 0)), 
		1 * kpc
	);
	
	double minStep = 0.1 * kpc;
	propa.setMinimumStep(minStep);
	propa.setTolerance(0.42);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step
	EXPECT_DOUBLE_EQ(5 * propa.getMinimumTimeStep(), c.getNextStep());  // acceleration by factor 5
}


#ifndef CRPROPA_TESTS_SKIP_EXCEPTIONS
TEST(testPropagationBP, exceptionsTimeStep) {
	// minStep should be smaller than maxStep
	EXPECT_THROW(PropagationBP propa(0.42, 10 , 0, new UniformMagneticField(Vector3d(0, 0, 1 * nG))), std::runtime_error);
	// Too large tolerance: tolerance should be between 0 and 1
	EXPECT_THROW(PropagationBP propa(42., 10 * kiloyear , 20 * kiloyear, new UniformMagneticField(Vector3d(0, 0, 1 * nG))), std::runtime_error);

	PropagationBP propa(1 * Megayear, new UniformMagneticField(Vector3d(0, 0, 1 * nG)));

	// set maximum step, so that it can be tested what happens if a larger minStep is set.
	propa.setMaximumTimeStep(1 * Megayear);

	// this tests _that_ the expected exception is thrown
	EXPECT_THROW(propa.setTolerance(2.), std::runtime_error);
	EXPECT_THROW(propa.setMinimumTimeStep(-1.), std::runtime_error);
	EXPECT_THROW(propa.setMinimumTimeStep(2 * Megayear), std::runtime_error);

	// set minimum step, so that it can be tested what happens if a smaller maxStep is set.
	propa.setMinimumTimeStep(0.5 * Megayear);

	EXPECT_THROW(propa.setMaximumTimeStep(0.1 * Megayear), std::runtime_error);
}


TEST(testPropagationBP, exceptions) {
	// minStep should be smaller than maxStep
	EXPECT_THROW(PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 0.42, 10 , 0), std::runtime_error);
	// Too large tolerance: tolerance should be between 0 and 1
	EXPECT_THROW(PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)), 42., 10 * kpc , 20 * kpc), std::runtime_error);

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
#endif


TEST(testPropagationBP, constructorTimeStep) {
	// Test construction and parameters
	// default values:
	ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(0, 0, 1 * nG));
	double minStep = 1.;
	double maxStep = 100.;
	double tolerance = 0.01;

	{  // check default constructor (only supports length steps)
		PropagationBP propa;

		EXPECT_FALSE(propa.getField().valid());
		EXPECT_EQ(1 * kiloparsec / c_light, propa.getMinimumTimeStep());
		EXPECT_EQ(1 * kiloparsec / c_light, propa.getMaximumTimeStep());
	}

	{  // test fixed stepsize constructor
		PropagationBP propa(minStep, bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(minStep, propa.getMaximumTimeStep());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test B-Field only constructor
		PropagationBP propa(tolerance, minStep, maxStep, bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(maxStep, propa.getMaximumTimeStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test fixed step constructor with B-Field and E-Field
		PropagationBP propa(minStep, bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(minStep, propa.getMaximumTimeStep());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test B-Field + E-Field constructor
		PropagationBP propa(tolerance, minStep, maxStep, bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(maxStep, propa.getMaximumTimeStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{	// test update of parameter
		PropagationBP propa;

		// Update parameters
		double minStep = 10.;
		double maxStep = 10.;
		double tolerance = 0.0001;
		ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(10 * nG, 0, 1 * nG));
		
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		propa.setTolerance(tolerance);
		propa.setField(bField);

		EXPECT_EQ(minStep, propa.getMinimumTimeStep());
		EXPECT_EQ(maxStep, propa.getMaximumTimeStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}
}


TEST(testPropagationBP, constructor) {
	// Test construction and parameters
	// default values:
	ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(0, 0, 1 * nG));
	double minStep = 1.;
	double maxStep = 100.;
	double tolerance = 0.01;

	{  // check default constructor
		PropagationBP propa;

		EXPECT_FALSE(propa.getField().valid());
		EXPECT_EQ(1 * kiloparsec, propa.getMinimumStep());
		EXPECT_EQ(1 * kiloparsec, propa.getMaximumStep());
	}

	{  // test fixed stepsize constructor
		PropagationBP propa(bField, minStep);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(minStep, propa.getMaximumStep());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test B-Field only constructor
		PropagationBP propa(bField, tolerance, minStep, maxStep);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(maxStep, propa.getMaximumStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test fixed step constructor with B-Field and E-Field
		PropagationBP propa(bField, minStep);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(minStep, propa.getMaximumStep());
		EXPECT_EQ(bField, propa.getField());
	}

	{  // test B-Field + E-Field constructor
		PropagationBP propa(bField, tolerance, minStep, maxStep);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(maxStep, propa.getMaximumStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}

	{	// test update of parameter
		PropagationBP propa;

		// Update parameters
		double minStep = 10.;
		double maxStep = 10.;
		double tolerance = 0.0001;
		ref_ptr<MagneticField> bField = new UniformMagneticField(Vector3d(10 * nG, 0, 1 * nG));
		
		propa.setMinimumStep(minStep);
		propa.setMaximumStep(maxStep);
		propa.setTolerance(tolerance);
		propa.setField(bField);

		EXPECT_EQ(minStep, propa.getMinimumStep());
		EXPECT_EQ(maxStep, propa.getMaximumStep());
		EXPECT_EQ(tolerance, propa.getTolerance());
		EXPECT_EQ(bField, propa.getField());
	}
}


// Test if the step size is reduced correctly if the error is too large with respect to the tolerance: r > 1
TEST(testPropagationBP, reduceTimeStep) {
	
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * TeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	
	{ // test reduction with only B-Field:
		Candidate c(p);
		PropagationBP propa(1. * kiloyear, new UniformMagneticField(Vector3d(0, 0, 100 * nG)));

		double minStep = 0.1 * kpc / c.getVelocity();
		double maxStep = 1 * Gpc / c.getVelocity();
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		// small tolerance leads to large values of r
		propa.setTolerance(1e-15);
		c.setNextStep(propa.getMaximumTimeStep());

		propa.process(&c);

		// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
	}

	{  // test with B-Field and E-Field:
		Candidate c(p);
		PropagationBP propa(
			1. * kiloyear, 
			new UniformMagneticField(Vector3d(0, 0, 100 * nG))
		);

		double minStep = 0.1 * kpc / c.getVelocity();
		double maxStep = 1 * Gpc / c.getVelocity();
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		// small tolerance leads to large values of r
		propa.setTolerance(1e-15);
		c.setNextStep(propa.getMaximumTimeStep());

		propa.process(&c);

		// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
	}

	{  // test with low energy:
		ParticleState p;
		p.setId(nucleusId(1, 1));
		p.setEnergy(1 * eV);
		p.setPosition(Vector3d(0, 0, 0));
		p.setDirection(Vector3d(0, 1, 0));
		
		Candidate c(p);
		PropagationBP propa(
			1. * kiloyear, 
			new UniformMagneticField(Vector3d(0, 0, 100 * milli * tesla))
		);

		double minStep = 0.1 * meter / p.getVelocity().getR();
		double maxStep = 1. * kilometer / p.getVelocity().getR();
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		// small tolerance leads to large values of r
		propa.setTolerance(1e-15);
		c.setNextStep(propa.getMaximumTimeStep());

		propa.process(&c);

		// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
	}
}


// Test if the step size is reduced correctly if the error is too large with respect to the tolerance: r > 1
TEST(testPropagationBP, reduceStep) {
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * TeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	
	{ // test reduction with only B-Field:
		Candidate c(p);
		PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 100 * nG)));

		double minStep = 0.1 * kpc;
		double maxStep = 1 * Gpc;
		propa.setMinimumStep(minStep);
		propa.setMaximumStep(maxStep);
		// small tolerance leads to large values of r
		propa.setTolerance(1e-15);
		c.setNextStep(propa.getMaximumTimeStep());

		propa.process(&c);

		// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
	}

	{  // test with B-Field and E-Field:
		Candidate c(p);
		PropagationBP propa(
			new UniformMagneticField(Vector3d(0, 0, 100 * nG)),
			1. * kpc
		);

		double minStep = 0.1 * kpc;
		double maxStep = 1 * Gpc;
		propa.setMinimumStep(minStep);
		propa.setMaximumStep(maxStep);
		// small tolerance leads to large values of r
		propa.setTolerance(1e-15);
		c.setNextStep(propa.getMaximumTimeStep());

		propa.process(&c);

		// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
	}

	{  // test with low energy:
		ParticleState p;
		p.setId(nucleusId(1, 1));
		p.setEnergy(1 * eV);
		p.setPosition(Vector3d(0, 0, 0));
		p.setDirection(Vector3d(0, 1, 0));
		
		Candidate c(p);
		PropagationBP propa(
			new UniformMagneticField(Vector3d(0, 0, 100 * milli * tesla)),
			1. * kpc
		);

		double minStep = 0.1 * meter / p.getVelocity().getR();
		double maxStep = 1. * kilometer / p.getVelocity().getR();
		// since it is converted internally with c_light and we are very far off:
		propa.setMinimumStep(minStep*c_light);
		propa.setMaximumStep(maxStep*c_light);
		// small tolerance leads to large values of r
		propa.setTolerance(1e-15);
		c.setNextStep(propa.getMaximumTimeStep());

		propa.process(&c);

		// adaptive algorithm should propagate particle with minimum step size due to the low value for the tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getCurrentStep());  // perform minimum step because of large r due to small tolerance
		EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep(), c.getNextStep());  // stay at minimum step because of large r due to small tolerance
	}
}


// Test if the step size is increased correctly if the error is small with respect to the tolerance: r < 1
TEST(testPropagationBP, increaseTimeStep) {
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	
	{  // test for only B-Field
		Candidate c(p);
		PropagationBP propa(1 * kiloyear, new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
		
		double minStep = 0.001 * pc/c.getVelocity();
		double maxStep = 3.125 * pc/c.getVelocity();
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		// large tolerance leads to small values of r. Consequently, the step size can be increased.
		propa.setTolerance(0.9);

		// each step the step size can be increased by a factor of 5.
		for (int i = 1; i < 6; i++){
			propa.process(&c);
			EXPECT_DOUBLE_EQ(minStep*pow(5, i), c.getNextStep());
		}
		// after 5 steps the maxStep is reached. The current step is, however, less.
		EXPECT_DOUBLE_EQ(maxStep/5., c.getCurrentStep());
		EXPECT_DOUBLE_EQ(maxStep, c.getNextStep());
	}

	{  // test for E-Field + B-Field
		Candidate c(p);
		PropagationBP propa(
			1 * kiloyear, 
			new UniformMagneticField(Vector3d(0, 0, 1 * nG))
		);
		
		double minStep = 0.001 * pc/c.getVelocity();
		double maxStep = 3.125 * pc/c.getVelocity();
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		// large tolerance leads to small values of r. Consequently, the step size can be increased.
		propa.setTolerance(0.9);

		// each step the step size can be increased by a factor of 5.
		for (int i = 1; i < 6; i++){
			propa.process(&c);
			EXPECT_DOUBLE_EQ(minStep*pow(5, i), c.getNextStep());
		}
		// after 5 steps the maxStep is reached. The current step is, however, less.
		EXPECT_DOUBLE_EQ(maxStep/5., c.getCurrentStep());
		EXPECT_DOUBLE_EQ(maxStep, c.getNextStep());
	}

	{  // test with low energy
		ParticleState p;
		p.setId(nucleusId(1, 1));
		p.setEnergy(1 * eV);
		p.setPosition(Vector3d(0, 0, 0));
		p.setDirection(Vector3d(0, 1, 0));
		Candidate c(p);
		
		PropagationBP propa(
			1 * kiloyear, 
			new UniformMagneticField(Vector3d(0, 0, 1 * nano*tesla))
		);
		
		double minStep = 0.001 * cm/c.getVelocity();
		double maxStep = 3.125 * cm/c.getVelocity();
		propa.setMinimumTimeStep(minStep);
		propa.setMaximumTimeStep(maxStep);
		// large tolerance leads to small values of r. Consequently, the step size can be increased.
		propa.setTolerance(0.9);

		// each step the step size can be increased by a factor of 5.
		for (int i = 1; i < 6; i++){
			propa.process(&c);
			EXPECT_DOUBLE_EQ(minStep*pow(5, i), c.getNextStep());
		}
		// after 5 steps the maxStep is reached. The current step is, however, less.
		EXPECT_DOUBLE_EQ(maxStep/5., c.getCurrentStep());
		EXPECT_DOUBLE_EQ(maxStep, c.getNextStep());
	}
}


// Test if the step size is increased correctly if the error is small with respect to the tolerance: r < 1
TEST(testPropagationBP, increaseStep) {
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	
	{  // test reduction with only B-Field:
		Candidate c(p);
		PropagationBP propa(
			new UniformMagneticField(Vector3d(0, 0, 1 * nG)),
			1 * kiloparsec
		);
		
		double minStep = 0.001 * pc;
		double maxStep = 3.125 * pc;
		propa.setMinimumStep(minStep);
		propa.setMaximumStep(maxStep);
		// large tolerance leads to small values of r. Consequently, the step size can be increased.
		propa.setTolerance(0.9);

		// each step the step size can be increased by a factor of 5.
		for (int i = 1; i < 6; i++){
			propa.process(&c);
			EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep()*pow(5, i), c.getNextStep());
		}
		// after 5 steps the maxStep is reached. The current step is, however, less.
		EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep()/5., c.getCurrentStep());
		EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep(), c.getNextStep());
	}

	{  // test with B-Field and E-Field:
		Candidate c(p);
		PropagationBP propa(
			new UniformMagneticField(Vector3d(0, 0, 1 * nG)),
			1 * kiloparsec
		);
		
		double minStep = 0.001 * pc;
		double maxStep = 3.125 * pc;
		propa.setMinimumStep(minStep);
		propa.setMaximumStep(maxStep);
		// large tolerance leads to small values of r. Consequently, the step size can be increased.
		propa.setTolerance(0.9);

		// each step the step size can be increased by a factor of 5.
		for (int i = 1; i < 6; i++){
			propa.process(&c);
			EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep()*pow(5, i), c.getNextStep());
		}
		// after 5 steps the maxStep is reached. The current step is, however, less.
		EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep()/5., c.getCurrentStep());
		EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep(), c.getNextStep());
	}

	{  // test with low energy
		ParticleState p;
		p.setId(nucleusId(1, 1));
		p.setEnergy(1 * eV);
		p.setPosition(Vector3d(0, 0, 0));
		p.setDirection(Vector3d(0, 1, 0));
		Candidate c(p);
		
		PropagationBP propa(
			1 * kiloyear, 
			new UniformMagneticField(Vector3d(0, 0, 1 * nano*tesla))
		);
		
		double minStep = 0.001 * cm/c.getVelocity();
		double maxStep = 3.125 * cm/c.getVelocity();
		// since it is converted internally with c_light and we are very far off:
		propa.setMinimumStep(minStep*c_light);
		propa.setMaximumStep(maxStep*c_light);
		// large tolerance leads to small values of r. Consequently, the step size can be increased.
		propa.setTolerance(0.9);

		// each step the step size can be increased by a factor of 5.
		for (int i = 1; i < 6; i++){
			propa.process(&c);
			EXPECT_DOUBLE_EQ(propa.getMinimumTimeStep()*pow(5, i), c.getNextStep());
		}
		// after 5 steps the maxStep is reached. The current step is, however, less.
		EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep()/5., c.getCurrentStep());
		EXPECT_DOUBLE_EQ(propa.getMaximumTimeStep(), c.getNextStep());
	}
}


TEST(testPropagationBP, protonTimeStep) {
	
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));
	
	PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
	Candidate c(p);
	c.setNextStep(0);

	double step = 0.001 * kpc/c.getVelocity();
	propa.setMinimumTimeStep(step);
	propa.setMaximumTimeStep(10*step);
	propa.setTolerance(0.00001);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(step, c.getCurrentStep());  // perform step
	EXPECT_DOUBLE_EQ(5 * step, c.getNextStep());  // acceleration by factor 5
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

	EXPECT_DOUBLE_EQ(step/c.getVelocity(), c.getCurrentStep());  // perform step
	EXPECT_DOUBLE_EQ(5 * step/c.getVelocity(), c.getNextStep());  // acceleration by factor 5
}


TEST(testPropagationBP, neutronTimeStep) {
	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	PropagationBP propa(
		new UniformMagneticField(Vector3d(0, 0, 1 * nG)),
		1 * kpc
	);
	
	double step = kpc/c.getVelocity();
	
	propa.setMaximumTimeStep(step);
	propa.setMinimumTimeStep(step);

	propa.process(&c);

	EXPECT_DOUBLE_EQ(step, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(step, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}


TEST(testPropagationBP, neutron) {
	ParticleState p;
	p.setId(nucleusId(1, 0));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 1, 0));

	Candidate c(p);
	PropagationBP propa(
		new UniformMagneticField(Vector3d(0, 0, 1 * nG)),
		1 * kpc
	);
	
	double step = kpc;
	
	propa.setMaximumStep(step);
	propa.setMinimumStep(step);

	propa.process(&c);

	// divide by c_light since propagator does internal conversion over 1/c_light
	EXPECT_DOUBLE_EQ(step/c_light, c.getCurrentStep());
	EXPECT_DOUBLE_EQ(step/c_light, c.getNextStep());
	EXPECT_EQ(Vector3d(0, 1 * kpc, 0), c.current.getPosition());
	EXPECT_EQ(Vector3d(0, 1, 0), c.current.getDirection());
}


// Test the numerical results for parallel magnetic field lines along the z-axis
TEST(testPropagationBP, gyrationTimeStep) {
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(100 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(1, 1, 1));
	
	Candidate c(p);
	PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
	
	double step = 10. * Mpc/c.getVelocity();  // gyroradius is 108.1 Mpc
	propa.setMaximumTimeStep(step);
	propa.setMinimumTimeStep(step);
	c.setNextStep(0);

	propa.process(&c);

	double dirX = c.current.getDirection().x;
	double dirY = c.current.getDirection().y;
	double dirZ = c.current.getDirection().z;
	double posZ = c.current.getPosition().z;

	// Test if the analytical solution is achieved for the components of the momentum with the Boris push as expected in
	// the background magnetic field.
	EXPECT_DOUBLE_EQ(2 / 3., dirX * dirX + dirY * dirY);  // constant momentum in the perpendicular plane to background magnetic field field
	EXPECT_DOUBLE_EQ(1 / 3., dirZ * dirZ);  // constant momentum parallel to the background magnetic field
	EXPECT_DOUBLE_EQ( step * step * c.getVelocity() * c.getVelocity() / 3., posZ * posZ);  // constant velocity parallel to the background magnetic field

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
	EXPECT_DOUBLE_EQ(100 * step * step * c.getVelocity() * c.getVelocity() / 3., posZ * posZ);  // constant velocity parallel to the background magnetic field
}


// Test the numerical results for parallel magnetic field lines along the z-axis
TEST(testPropagationBP, gyration) {
	PropagationBP propa(new UniformMagneticField(Vector3d(0, 0, 1 * nG)));
	
	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(1000 * EeV);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(1, 1, 1));
	Candidate c(p);
	
	double step = 10. * Mpc;  // gyroradius is 108.1 Mpc
	propa.setMaximumStep(step);
	propa.setMinimumStep(step);
	c.setNextStep(0);

	propa.process(&c);

	double dirX = c.current.getDirection().x;
	double dirY = c.current.getDirection().y;
	double dirZ = c.current.getDirection().z;
	double posZ = c.current.getPosition().z;

	// Test if the analytical solution is achieved for the components of the momentum with the Boris push as expected in
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
	// it is only expected to be near 1/3. since we hand over stepsize for a particle with velocity only approx. c_light
	EXPECT_DOUBLE_EQ(1 / 3., dirZ * dirZ);  // constant momentum parallel to the background magnetic field
	EXPECT_DOUBLE_EQ(100 * step * step / 3., posZ * posZ);  // constant velocity parallel to the background magnetic field
}


// Test the that the optimization for fixed step sizes works
TEST(testPropagationBP, fixedTimeStepOptimization) {
	// particle 1 with fixed step sizes
	double fixed_step = kiloyear;
	PropagationBP propa1(fixed_step, new PlaneWaveTurbulence(TurbulenceSpectrum(gauss, pc, 100*pc), 10, 1));
	ParticleState p1;
	p1.setId(nucleusId(1, 1));
	p1.setEnergy(100 * EeV);
	p1.setPosition(Vector3d(0, 0, 0));
	p1.setDirection(Vector3d(1, 1, 1));
	Candidate c1(p1);
	c1.setNextStep(0);
	// Nine new steps to have finally propagated the particle ten times
	for (int i = 0; i < 9; i++){
		propa1.process(&c1);
	}

	// particle 2 with different min and max steps. The tolerance is chosen such that particle 2 will be
	// propagated with the same step as particle 1, however not using the optimization for fixed step sizes
	double tolerance = 0.001;
	PropagationBP propa2(tolerance, fixed_step, 1.1*fixed_step, new PlaneWaveTurbulence(TurbulenceSpectrum(gauss, pc, 100*pc), 10, 1));
	ParticleState p2;
	p2.setId(nucleusId(1, 1));
	p2.setEnergy(100 * EeV);
	p2.setPosition(Vector3d(0, 0, 0));
	p2.setDirection(Vector3d(1, 1, 1));
	Candidate c2(p2);
	c1.setNextStep(0);
	// Nine new steps to have finally propagated the particle ten times
	for (int i = 0; i < 9; i++){
		propa2.process(&c2);
	}

	EXPECT_DOUBLE_EQ(c1.current.getDirection().x, c2.current.getDirection().x);
	EXPECT_DOUBLE_EQ(c1.current.getDirection().y, c2.current.getDirection().y);
	EXPECT_DOUBLE_EQ(c1.current.getDirection().z, c2.current.getDirection().z);
	EXPECT_DOUBLE_EQ(c1.current.getPosition().x, c2.current.getPosition().x);
	EXPECT_DOUBLE_EQ(c1.current.getPosition().y, c2.current.getPosition().y);
	EXPECT_DOUBLE_EQ(c1.current.getPosition().z, c2.current.getPosition().z);
}


// Test the that the optimization for fixed step sizes works
TEST(testPropagationBP, fixedStepOptimization) {
	// particle 1 with fixed step sizes
	double fixed_step = pc;
	PropagationBP propa1(new PlaneWaveTurbulence(TurbulenceSpectrum(gauss, pc, 100*pc), 10, 1), fixed_step);
	ParticleState p1;
	p1.setId(nucleusId(1, 1));
	p1.setEnergy(100 * EeV);
	p1.setPosition(Vector3d(0, 0, 0));
	p1.setDirection(Vector3d(1, 1, 1));
	Candidate c1(p1);
	c1.setNextStep(0);
	// Nine new steps to have finally propagated the particle ten times
	for (int i = 0; i < 9; i++){
		propa1.process(&c1);
	}

	// particle 2 with different min and max steps. The tolerance is chosen such that particle 2 will be
	// propagated with the same step as particle 1, however not using the optimization for fixed step sizes
	double tolerance = 0.001;
	PropagationBP propa2(new PlaneWaveTurbulence(TurbulenceSpectrum(gauss, pc, 100*pc), 10, 1), tolerance, fixed_step, 1.1*fixed_step);
	ParticleState p2;
	p2.setId(nucleusId(1, 1));
	p2.setEnergy(100 * EeV);
	p2.setPosition(Vector3d(0, 0, 0));
	p2.setDirection(Vector3d(1, 1, 1));
	Candidate c2(p2);
	c1.setNextStep(0);
	// Nine new steps to have finally propagated the particle ten times
	for (int i = 0; i < 9; i++){
		propa2.process(&c2);
	}

	EXPECT_DOUBLE_EQ(c1.current.getDirection().x, c2.current.getDirection().x);
	EXPECT_DOUBLE_EQ(c1.current.getDirection().y, c2.current.getDirection().y);
	EXPECT_DOUBLE_EQ(c1.current.getDirection().z, c2.current.getDirection().z);
	EXPECT_DOUBLE_EQ(c1.current.getPosition().x, c2.current.getPosition().x);
	EXPECT_DOUBLE_EQ(c1.current.getPosition().y, c2.current.getPosition().y);
	EXPECT_DOUBLE_EQ(c1.current.getPosition().z, c2.current.getPosition().z);
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

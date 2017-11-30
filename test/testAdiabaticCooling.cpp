#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/advectionField/AdvectionField.h"
#include "crpropa/module/AdiabaticCooling.h"
#include "gtest/gtest.h"

//#include <fstream>

namespace crpropa {

// AdiabaticCooling ---------------------------------------------------------------

TEST (AdiabaticCooling, UniformField) {
	// Test in a uniform advection Field

	AdiabaticCooling AC(new UniformAdvectionField(Vector3d(1,0,0)));
	Candidate c(nucleusId(1,1), 1e13*eV);
	c.setCurrentStep(10*kpc);
	c.setNextStep(10*kpc);
	double E = c.current.getEnergy();
	AC.process(&c);

	// Energy is expected to be conserved
	EXPECT_DOUBLE_EQ(c.current.getEnergy(), E);
	EXPECT_DOUBLE_EQ(c.getNextStep(), 10*kpc);

	double limit = 0.2;
	AdiabaticCooling AC2(new UniformAdvectionField(Vector3d(1,0,0)), limit);
	
	EXPECT_DOUBLE_EQ(AC2.getLimit(), limit);

	//

}

TEST (AdiabaticCooling, ConstantSphericalField) {
	// Constant velocity vector
	
	AdiabaticCooling AC(new ConstantSphericalAdvectionField(Vector3d(0,0,0), 1));
	Candidate c(nucleusId(1,1), 10);
	c.current.setPosition(Vector3d(1,0,0));
	c.setCurrentStep(c_light);
	c.setNextStep(c_light);
	double E = c.current.getEnergy();
	AC.process(&c);

	// Check energy loss and step limitation
	EXPECT_DOUBLE_EQ(c.current.getEnergy(), E/3.);
	EXPECT_DOUBLE_EQ(c.getNextStep(), 0.15*c_light);

}


} // namespace crpropa

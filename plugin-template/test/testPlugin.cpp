#include "CRPropa.h"
#include "myPlugin.h"

#include "gtest/gtest.h"



// With this expression we do not need to write `crpropa::` for everything out of the crpropa namespace
using namespace crpropa;

// create a test over TEST(ModuleToTest, WhatToTest){}
TEST(MyModule, SimpleTest) {

	ModuleList sim;
	sim.add(new SimplePropagation(1 * pc, 1 * pc));
	sim.add(new MaximumTrajectoryLength(1000 * pc));
	sim.add(new myPlugin::MyModule());

	// define source to test this module and get candidate with correct flag out of there
	Source source;

	// `AddMyProperty` is added in plugin.h but under the crpropa namespace:
	source.add(new AddMyProperty());
	
	// check if the candidate has the correct property set by the source
	ref_ptr<Candidate> cand = source.getCandidate();
	EXPECT_EQ(cand->getProperty("counter"), uint32_t(0));

	sim.run(cand);

	// in a test you can test your values like this:
	EXPECT_EQ(cand->getProperty("counter"), uint32_t(1000));
	// but there are many more such EXPECT expressions
}

int main(int argc, char **argv) {
	// this runs all above defined tests
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
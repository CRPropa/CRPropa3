/** Unit tests for output modules */

#include "crpropa/module/OutputTXT.h"
#include "crpropa/module/OutputCRPropa2.h"
#include "crpropa/module/OutputROOT.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(ConditionalOutput, removeProperty) {
	ConditionalOutput output("", "Detected");
	Candidate c;
	c.setProperty("Dected", "");
	output.process(&c);
	EXPECT_FALSE(c.hasProperty("Detected"));
}

TEST(EventOutput1D, removeProperty) {
	EventOutput1D output("");
	Candidate c;
	c.setProperty("Dected", "");
	output.process(&c);
	EXPECT_FALSE(c.hasProperty("Detected"));
}

TEST(CRPropa2EventOutput3D, removeProperty) {
	CRPropa2EventOutput3D output("");
	Candidate c;
	c.setProperty("Dected", "");
	output.process(&c);
	EXPECT_FALSE(c.hasProperty("Detected"));
}

TEST(CRPropa2EventOutput1D, removeProperty) {
	CRPropa2EventOutput1D output("");
	Candidate c;
	c.setProperty("Dected", "");
	output.process(&c);
	EXPECT_FALSE(c.hasProperty("Detected"));
}

#ifdef CRPROPA_HAVE_ROOT
TEST(CRPropa2ROOTEventOutput3D, removeProperty) {
	CRPropa2ROOTEventOutput3D output("");
	Candidate c;
	c.setProperty("Dected", "");
	output.process(&c);
	EXPECT_FALSE(c.hasProperty("Detected"));
}

TEST(CRPropa2ROOTEventOutput1D, removeProperty) {
	CRPropa2ROOTEventOutput1D output("");
	Candidate c;
	c.setProperty("Dected", "");
	output.process(&c);
	EXPECT_FALSE(c.hasProperty("Detected"));
}
#endif

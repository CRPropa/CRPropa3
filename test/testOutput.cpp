/** Unit tests for Output modules of CRPropa
    Output
    TextOutput
    ParticleCollector
 */

#include "crpropa/ParticleID.h"
#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/Version.h"
#include "crpropa/module/Output.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/module/ParticleCollector.h"

#include <string>
#include "gtest/gtest.h"
#include <iostream>

namespace crpropa {

//-- Output

TEST(Output, getCount) {
	Candidate c;
	Output output;
	for (int it=0; it<5; ++it, output.process(&c));

	EXPECT_EQ(output.getCount(), 5);
}

//-- TextOutput

TEST(TextOutput, printHeader_Trajectory1D) {
	Candidate c;
	TextOutput output(Output::Trajectory1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tID\tE\tX");
}

TEST(TextOutput, printHeader_Event1D) {
	Candidate c;
	TextOutput output(Output::Event1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tD\tID\tE\tID0\tE0");
}

TEST(TextOutput, printHeader_Trajectory3D) {
	Candidate c;
	TextOutput output(Output::Trajectory3D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tD\tID\tE\tX\tY\tZ\tPx\tPy\tPz");
}

TEST(TextOutput, printHeader_Event3D) {
	Candidate c;
	TextOutput output(Output::Event3D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tD\tID\tE\tX\tY\tZ\tPx\tPy\tPz\tID0\tE0\tX0\tY0\tZ0\tP0x\tP0y\tP0z");
}

TEST(TextOutput, printHeader_Custom) {
	Candidate c;
	TextOutput output(Output::Event1D);

	output.enable(Output::SerialNumberColumn);
	output.disable(Output::TrajectoryLengthColumn);
	output.set(Output::RedshiftColumn, false);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tSN\tID\tE\tSN0\tID0\tE0\tSN1");
}

TEST(TextOutput, printProperty) {
	Candidate c;
	TextOutput output(Output::Event1D);
	output.disableAll();
	output.enableProperty("foo", 2.0, "Bar");

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	// name in first line of header
	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tfoo");
}

TEST(TextOutput, printHeader_Version) {
	Candidate c;
	TextOutput output(Output::Event1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	// length of the prefix is 19 chars
	size_t version_pos = captured.find("# CRPropa version: ") + 19;

	EXPECT_EQ(
		captured.substr(
			version_pos,
		      	captured.find("\n", version_pos) - version_pos
			),
	         g_GIT_DESC);
}

//-- ParticleCollector

TEST(ParticleCollector, getCount) {
	ref_ptr<Candidate> c = new Candidate();
	ParticleCollector output;

	for (int it=0; it<5; ++it, output.process(c));

	EXPECT_EQ(output.getCount(), 5);
}

TEST(ParticleCollector, fetchItem) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1,1), 1*EeV);
	ParticleCollector output;

	output.process(c);

	EXPECT_EQ(output[0], c);
}

TEST(ParticleCollector, reprocess) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1,1), 1*EeV);
	ParticleCollector collector;
	ParticleCollector output;

	collector.process(c);
	collector.reprocess(&output);

	EXPECT_EQ(output[0], c);
}

TEST(ParticleCollector, dumpload) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1,1), 1.234*EeV);
	c->current.setPosition(Vector3d(1,2,3));
	c->current.setDirection(Vector3d(-1,-1,-1));
	c->setTrajectoryLength(1*Mpc);
	c->setRedshift(2);

	ParticleCollector input;
	ParticleCollector output;

	for(int i=0; i<=10; ++i){
		input.process(c);
	}

	input.dump("ParticleCollector_DumpTest.txt");
	output.load("ParticleCollector_DumpTest.txt");

	EXPECT_EQ(input.getCount(), output.getCount());
	EXPECT_EQ(output[0]->current.getEnergy(), c->current.getEnergy());
	EXPECT_EQ(output[1]->getTrajectoryLength(), c->getTrajectoryLength());
	EXPECT_EQ(output[2]->current.getId(), c->current.getId());
	EXPECT_EQ(output[3]->getRedshift(), c->getRedshift());
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

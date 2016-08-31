/** Unit tests for Output modules of CRPropa
	Output
  	TextOutput
 */

#include "crpropa/ParticleID.h"
#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/module/Output.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/module/ParticleContainerOutput.h"

#include <string>
#include "gtest/gtest.h"

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
	          "#\tD\tID\tE\tX\tY\tZ\tPx\tPy\tPz\tID0\tE0\tX0\tY0\tZ0");
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

//-- ParticleContainerOutput

TEST(ParticleContainerOutput, getCount) {
	ref_ptr<Candidate> c = new Candidate();
	ParticleContainerOutput output;

	for (int it=0; it<5; ++it, output.process(c));
	
	EXPECT_EQ(output.getCount(), 5);
}

TEST(ParticleContainerOutput, fetchItem) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1,1), 1*EeV);
	ParticleContainerOutput output;

	output.process(c);
	
	EXPECT_EQ(output[0]->current.getEnergy(), 1*EeV);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

/** Unit tests for Output modules of CRPropa
    Output
    TextOutput
    ParticleCollector
 */

#include "CRPropa.h"

#include "gtest/gtest.h"
#include <iostream>
#include <string>


#ifdef CRPROPA_HAVE_HDF5
#include <hdf5.h>
#endif

// compare two arrays (intead of using Google Mock)
// https://stackoverflow.com/a/10062016/6819103
template <typename T, size_t size>
::testing::AssertionResult ArraysMatch(const T (&expected)[size],
                                       const T (&actual)[size]) {
	for (size_t i(0); i < size; ++i) {
		if (expected[i] != actual[i]) {
			return ::testing::AssertionFailure()
			       << "array[" << i << "] (" << actual[i] << ") != expected["
			       << i << "] (" << expected[i] << ")";
		}
	}

	return ::testing::AssertionSuccess();
}

namespace crpropa {

//-- Output

TEST(Output, size) {
	Candidate c;
	Output output;
	for (int it = 0; it < 5; ++it, output.process(&c));

	EXPECT_EQ(output.size(), 5);
}

//-- TextOutput

TEST(TextOutput, printHeader_Trajectory1D) {
	Candidate c;
	TextOutput output(Output::Trajectory1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")), "#\tID\tE\tX");
}

TEST(TextOutput, printHeader_Event1D) {
	Candidate c;
	TextOutput output(Output::Event1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	EXPECT_EQ(captured.substr(0, captured.find("\n")), "#\tD\tID\tE\tID0\tE0\ttag");
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

	EXPECT_EQ(
	    captured.substr(0, captured.find("\n")),
	    "#\tD\tID\tE\tX\tY\tZ\tPx\tPy\tPz\tID0\tE0\tX0\tY0\tZ0\tP0x\tP0y\tP0z\ttag");
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
	          "#\tSN\tID\tE\tSN0\tID0\tE0\tSN1\ttag");
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
	EXPECT_EQ(captured.substr(0, captured.find("\n")), "#\tfoo");
}

TEST(TextOutput, printHeader_Version) {
	Candidate c;
	TextOutput output(Output::Event1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	// length of the prefix is 19 chars
	size_t version_pos = captured.find("# CRPropa version: ") + 19;

	EXPECT_EQ(captured.substr(version_pos,
	                          captured.find("\n", version_pos) - version_pos),
	          g_GIT_DESC);
}

TEST(TextOutput, failOnIllegalOutputFile) {
	EXPECT_THROW(
	    TextOutput output("THIS_FOLDER_MUST_NOT_EXISTS_12345+/FILE.txt"),
	    std::runtime_error);
}

#ifdef CRPROPA_HAVE_HDF5
TEST(HDF5Output, failOnIllegalOutputFile) {
	HDF5Output out;
	// disable default error output of HDF5
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	EXPECT_THROW(out.open("THIS_FOLDER_MUST_NOT_EXISTS_12345+/FILE.h5"),
	             std::runtime_error);
}
#endif

//-- ParticleCollector

TEST(ParticleCollector, size) {
	ref_ptr<Candidate> c = new Candidate();
	ParticleCollector output;

	for (int it = 0; it < 5; ++it, output.process(c))
		;

	EXPECT_EQ(output.size(), 5);
}

TEST(ParticleCollector, fetchItem) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1, 1), 1 * EeV);
	ParticleCollector output;

	output.process(c);

	EXPECT_EQ(output[0], c);
}

TEST(ParticleCollector, reprocess) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1, 1), 1 * EeV);
	ParticleCollector collector;
	ParticleCollector output;

	collector.process(c);
	collector.reprocess(&output);

	EXPECT_EQ(output[0], c);
}

TEST(ParticleCollector, dumpload) {
	ref_ptr<Candidate> c = new Candidate(nucleusId(1, 1), 1.234 * EeV);
	c->current.setPosition(Vector3d(1, 2, 3));
	c->current.setDirection(Vector3d(-1, -1, -1));
	c->setTrajectoryLength(1 * Mpc);
	c->setRedshift(2);

	ParticleCollector input;
	ParticleCollector output;

	for (int i = 0; i <= 10; ++i) {
		input.process(c);
	}

	// Well, it would be nicer if we don't need to receate any file
	input.dump("ParticleCollector_DumpTest.txt");
	output.load("ParticleCollector_DumpTest.txt");

	EXPECT_EQ(input.size(), output.size());
	EXPECT_EQ(output[0]->current.getEnergy(), c->current.getEnergy());
	EXPECT_EQ(output[1]->getTrajectoryLength(), c->getTrajectoryLength());
	EXPECT_EQ(output[2]->current.getId(), c->current.getId());
	EXPECT_EQ(output[3]->getRedshift(), c->getRedshift());
}

// Just test if the trajectory is on a line for rectilinear propagation
TEST(ParticleCollector, getTrajectory) {
	int pos_x[10];
	int pos_x_expected[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

	ParticleState p;
	p.setPosition(Vector3d(10, 0, 0));
	p.setDirection(Vector3d(-1, 0, 0));
	ref_ptr<Candidate> c = new Candidate(p);

	ref_ptr<ParticleCollector> output = new ParticleCollector();
	ref_ptr<ParticleCollector> trajectory = new ParticleCollector();
	trajectory->setClone(true);

	ref_ptr<ModuleList> sim = new ModuleList();
	sim->add(new SimplePropagation(1, 1));

	ref_ptr<Observer> obs = new Observer();
	obs->add(new ObserverPoint());
	obs->onDetection(output);
	sim->add(obs);

	sim->run(c);

	output->getTrajectory(sim, 0, trajectory);

	Vector3d pos;
	int i = 0;

	for (ParticleCollector::iterator itr = trajectory->begin();
	     itr != trajectory->end(); ++itr) {
		pos = (*(itr->get())).current.getPosition();
		pos_x[i] = pos.getX();
		++i;
	}

	EXPECT_TRUE(ArraysMatch(pos_x_expected, pos_x));
}

TEST(ParticleCollector, runModuleList) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	modules.add(new MaximumTrajectoryLength(1 * Mpc));

	ParticleState p;
	p.setPosition(Vector3d(10, 0, 0));
	p.setDirection(Vector3d(-1, 0, 0));
	ref_ptr<Candidate> c = new Candidate(p);

	ref_ptr<ParticleCollector> collector = new ParticleCollector();

	collector->process(c);

	modules.setShowProgress(false);
	auto candidates = collector->getContainer();
	modules.run(&candidates);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

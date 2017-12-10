#include "crpropa/ModuleList.h"
#include "crpropa/Source.h"
#include "crpropa/ParticleID.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/BreakCondition.h"

#include "gtest/gtest.h"

namespace crpropa {

TEST(ModuleList, process) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	ParticleState initial;
	ref_ptr<Candidate> candidate = new Candidate(initial);
	modules.process(candidate);
}

TEST(ModuleList, getModule) {
	ModuleList modules;
	ref_ptr<SimplePropagation> prop = new SimplePropagation();
	modules.add(prop);
	EXPECT_TRUE(modules[0] == prop);
}

TEST(ModuleList, removeModule) {
	ModuleList modules;
	ref_ptr<SimplePropagation> prop = new SimplePropagation();
	modules.add(prop);
	modules.remove(0);
	EXPECT_EQ(modules.size(), 0);
}

TEST(ModuleList, runCandidateList) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	modules.add(new MaximumTrajectoryLength(1 * Mpc));
	ParticleState initial;
	ref_ptr<Candidate> candidate = new Candidate(initial);
	modules.run(candidate);
	EXPECT_DOUBLE_EQ(1 * Mpc, candidate->getTrajectoryLength());
	EXPECT_TRUE(candidate->isActive() == false);
}

TEST(ModuleList, runSource) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	modules.add(new MaximumTrajectoryLength(1 * Mpc));
	Source source;
	source.add(new SourcePosition(Vector3d(10, 0, 0) * Mpc));
	source.add(new SourceIsotropicEmission());
	source.add(new SourcePowerLawSpectrum(5 * EeV, 100 * EeV, -2));
	source.add(new SourceParticleType(nucleusId(1, 1)));
	modules.setShowProgress(true);
	modules.run(&source, 100, false);
}

#if _OPENMP
#include <omp.h>
TEST(ModuleList, runOpenMP) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	modules.add(new MaximumTrajectoryLength(1 * Mpc));
	Source source;
	source.add(new SourcePosition(Vector3d(10, 0, 0) * Mpc));
	source.add(new SourceIsotropicEmission());
	source.add(new SourcePowerLawSpectrum(5 * EeV, 100 * EeV, -2));
	source.add(new SourceParticleType(nucleusId(1, 1)));
	omp_set_num_threads(2);
	modules.run(&source, 1000, false);
}
#endif

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

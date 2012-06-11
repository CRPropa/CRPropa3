#include "mpc/ModuleList.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/Source.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(ModuleList, process) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	ParticleState initial;
	ref_ptr<Candidate> candidate = new Candidate(initial);
	modules.process(candidate);
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
	CompositeSource source(Vector3d(0.), 10., 100., -1.);
	source.addToComposition(getNucleusId(56, 26), 1);
	modules.setShowProgress(true);
	modules.run(&source, 100, false);
}

#if _OPENMP
#include <omp.h>
TEST(ModuleList, runOpenMP) {
	ModuleList modules;
	modules.add(new SimplePropagation());
	modules.add(new MaximumTrajectoryLength(1 * Mpc));
	CompositeSource source(Vector3d(0.), 10., 100., -1.);
	source.addToComposition(getNucleusId(56, 26), 1);
	omp_set_num_threads(2);
	modules.run(&source, 1000, false);
}
#endif

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

#include "crpropa/Candidate.h"
#include "crpropa/module/PhotonDINT1D.h"

#include "gtest/gtest.h"

namespace crpropa {

TEST(testDINT, step) {
	PhotonDINT1D mod("dint1d_output.txt");
	Candidate c;
	c.current.setId(22);
	c.current.setPosition(Vector3d(0.0001 * Mpc, 0, 0));
	c.current.setEnergy(10 * EeV);
	c.setActive(true);
	mod.process(&c);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

#include "gtest/gtest.h"
#include "mpc/DeflectionCK.h"



//TEST(testDeflectionCK, noBfield) {
//	HomogeneousMagneticField field(Hep3Vector(0., 0., 0.));
//	DeflectionCK d(1e-3);
//
//	Particle particle;
//	particle.setChargeMassNumber(1,1);
//	particle.setStepMpc(0.001);
//	particle.setNextStepMpc(0.001);
//	particle.setPositionMpc(Hep3Vector(0.,0.,0.));
//	particle.setDirection(Hep3Vector(0.,1.,0.));
//	particle.setEnergyEeV(1);
//
//	for (int i = 0; i < 10; i++)
//		d.apply(particle, field);
//	EXPECT_EQUAL


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

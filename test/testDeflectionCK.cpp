#include "gtest/gtest.h"
#include "DeflectionCK.h"


TEST(testDeflectionCK, cashKarpCoefficients) {
	// Runge-Kutta coefficients have to add up to 1 in every row of the Butcher table
	double sum;

	sum = 0;
	for (size_t i=0; i<6; i++)
		sum += cash_karp_b[i];
	EXPECT_DOUBLE_EQ(1, sum);

	sum = 0;
	for (size_t i=0; i<6; i++)
		sum += cash_karp_bs[i];
	EXPECT_DOUBLE_EQ(1, sum);

	for (size_t j=0; j<6; j++) {
		sum = 0;
		for (size_t i=0; i<5; i++)
			sum += cash_karp_a[i+j*5];
		EXPECT_DOUBLE_EQ(cash_karp_c[j], sum);
	}

}

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

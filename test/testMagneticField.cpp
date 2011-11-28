#include "gtest/gtest.h"
#include "mpc/MagneticField.h"
#include "mpc/TurbulentMagneticField.h"
#include <iostream>


namespace mpc {

TEST(testTurbulentMagneticField, ) {
	TurbulentMagneticField B(Vector3(0.), 64, 1, 1, -11. / 3., 2., 8.);
	B.initialize();
	a = field.Bx.real.mean().__abs__() <= precision
	b = field.By.real.mean().__abs__() <= precision
	c = field.Bz.real.mean().__abs__() <= precision
	self.assertTrue(a and b and c)
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

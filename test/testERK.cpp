#include "gtest/gtest.h"
#include "ExplicitRungeKutta.h"
#include "PhasePoint.h"
#include <iostream>

const double GM = 6.674e-11*5.974e24;
const double vr = 3075;
const double radius = 42157e3;

class CentralPotential: public ExplicitRungeKutta<PhasePoint>::F {
public:

	PhasePoint operator()(double t, const PhasePoint &v) {
		Hep3Vector acceleration = - v.a.unit() * GM / v.a.mag2();
		return PhasePoint(v.b, acceleration);
	}
};

//TEST(erkTest, cashCarpCoefficients) {
//	for (size_t i=0; i<6, i++) {
//		cash_karp_a[i];
//	}
//}

TEST(erkTest, position) {
	PhasePoint p(Hep3Vector(radius,0,0), Hep3Vector(0,vr,0)), pOut, pErr;
	ExplicitRungeKutta<PhasePoint> erk;
	erk.loadCashKarp();
	CentralPotential f;
	for (size_t i = 0; i<24; i++) {
		pErr.a = Hep3Vector(0,0,0);
		erk.step(0, p, pOut, pErr, 1, f); // step of 1 second
		std::cout << pErr.a/1000 << ", " << pErr.b << std::endl;
		p = pOut;
	}
	std::cout << p.a/1000 << ", " << p.b << std::endl;
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#include "gtest/gtest.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/PhasePoint.h"
#include <iostream>


namespace mpc {

TEST(testERK, cashKarpCoefficients) {
	// Runge-Kutta coefficients have to add up to 1 in every row of the Butcher table
	double sum;

	sum = 0;
	for (size_t i = 0; i < 6; i++)
		sum += cash_karp_b[i];
	EXPECT_DOUBLE_EQ(1, sum);

	sum = 0;
	for (size_t i = 0; i < 6; i++)
		sum += cash_karp_bs[i];
	EXPECT_DOUBLE_EQ(1, sum);

	for (size_t j = 0; j < 6; j++) {
		sum = 0;
		for (size_t i = 0; i < 5; i++)
			sum += cash_karp_a[i + j * 5];
		EXPECT_DOUBLE_EQ(cash_karp_c[j], sum);
	}
}

class CentralPotential: public ExplicitRungeKutta<PhasePoint>::F {
public:
	PhasePoint operator()(double t, const PhasePoint &v) {
		double GM = 3.986004418e14; // gravitational constant * mass of earth
		Vector3 acceleration = -v.a.unit() * GM / v.a.mag2();
		return PhasePoint(v.b, acceleration);
	}
};

TEST(erkTest, geoStationaryOrbit) {
	double siderialDay = 86164.0954; // siderial day in seconds
	double GM = 3.986004418e14; // gravitational constant * mass of earth
	double f = 2*3.14159/siderialDay; // angular frequency
	double r = pow(GM/(f*f), 1./3.); // radius
	double v = r*f; // radial velocity
	ExplicitRungeKutta<PhasePoint> erk;
	erk.loadCashKarp();
	PhasePoint y(Vector3(r, 0, 0), Vector3(0, v, 0)), yOut, yErr;
	CentralPotential dydt;

	double step = 100; // time step in seconds
	for (size_t i = 0; i < int(siderialDay/step); i++) {
		erk.step(0, y, yOut, yErr, step, dydt); // step of 1 second
		y = yOut;
		std::cout << y.a.x() << "," << y.a.y() << "," << y.a.z() << std::endl;
//		std::cout << y.a << ", " << y.b << std::endl;
	}
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

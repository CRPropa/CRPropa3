#include "dint/gauleg.h"

#include "gtest/gtest.h"
#include <math.h>

namespace crpropa {

class GaussLegendreIntegrationTest : public testing::Test {
protected:
	size_t N;
	double *x;
	double *w;
	double x1, x2;
	double analytical, res;
	virtual void SetUp() {
		N = 16384;
		x = new double[N];
		w = new double[N];
	}

	virtual void TearDown() {
		delete[] x;
		delete[] w;
	}
};

TEST_F(GaussLegendreIntegrationTest, exp) {
	x1 = -1;
	x2 = 1;

	res = 0;
	Gauleg(x1, x2, x, w, N);
	for (int i =0; i < N; i++)
		res+= w[i] * exp(x[i]);

	analytical = exp(x2) - exp(x1);
	EXPECT_NEAR(res, analytical, 1E-10);
}

TEST_F(GaussLegendreIntegrationTest, x2) {
	x1 = -10;
	x2 = 1;

	res = 0;
	Gauleg(x1, x2, x, w, N);
	for (int i =0; i < N; i++)
		res+= w[i] * x[i] * x[i];

	analytical = 1./3 * (x2*x2*x2 - x1*x1*x1);
	EXPECT_NEAR(res, analytical, 1E-10);
}

TEST_F(GaussLegendreIntegrationTest, x3) {
	x1 = -4;
	x2 = 7;

	res = 0;
	Gauleg(x1, x2, x, w, N);
	for (int i =0; i < N; i++)
		res+= w[i] * x[i] * x[i] * x[i];

	analytical = 1./4 * (x2*x2*x2*x2 - x1*x1*x1*x1);
	EXPECT_NEAR(res, analytical, 1E-10);
}

TEST_F(GaussLegendreIntegrationTest, cosx) {
	x1 = -2;
	x2 = 3;

	res = 0;
	Gauleg(x1, x2, x, w, N);
	for (int i =0; i < N; i++)
		res+= w[i] * cos(x[i] * 1234);

	analytical = 1./1234 * (sin(x2 * 1234 ) - sin(x1 * 1234));
	EXPECT_NEAR(res, analytical, 1E-10);
}

TEST_F(GaussLegendreIntegrationTest, sinx) {
	x1 = -M_PI;
	x2 = M_PI;

	res = 0;
	Gauleg(x1, x2, x, w, N);
	for (int i =0; i < N; i++)
		res+= w[i] * sin(x[i]);

	analytical = 0;
	EXPECT_NEAR(res, analytical, 1E-10);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

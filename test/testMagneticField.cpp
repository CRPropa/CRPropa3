#include "gtest/gtest.h"
#include "mpc/magneticfield/MagneticField.h"
#include "mpc/magneticfield/turbulentMagneticField.h"
#include "mpc/Vector3.h"
#include "fftw3.h"
#include <iostream>

namespace mpc {

TEST(testTurbulentMagneticField, Bmean) {
	// <B> = 0

	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();
	size_t n = B.n;

	Vector3 b(0, 0, 0);
	for (unsigned int ix = 0; ix < n; ix++)
		for (unsigned int iy = 0; iy < n; iy++)
			for (unsigned int iz = 0; iz < n; iz++)
				b += B.field[ix][iy][iz];
	b /= n * n * n;

	double precision = 1e-12;
	EXPECT_NEAR(b.x(), 0, precision);
	EXPECT_NEAR(b.y(), 0, precision);
	EXPECT_NEAR(b.z(), 0, precision);
}

TEST(testTurbulentMagneticField, Brms) {
	// <B^2> = Brms^2

	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();
	size_t n = B.n;

	double brms = 0;
	for (unsigned int ix = 0; ix < n; ix++)
		for (unsigned int iy = 0; iy < n; iy++)
			for (unsigned int iz = 0; iz < n; iz++)
				brms += B.field[ix][iy][iz].mag2();
	brms = sqrt(brms / n / n / n);

	double precision = 1e-12;
	EXPECT_NEAR(brms, 1, precision);
}

TEST(testTurbulentMagneticField, Solenoidity) {
	// B(k)*k = 0

	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();
	size_t n = B.n;

	// transform back into configuration space
	fftw_complex *Bx, *By, *Bz;
	Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
	By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
	Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);

	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				int i = ix * n * n + iy * n + iz;
				Bx[i][0] = B.field[ix][iy][iz].x();
				Bx[i][1] = 0;
				By[i][0] = B.field[ix][iy][iz].y();
				By[i][1] = 0;
				Bz[i][0] = B.field[ix][iy][iz].z();
				Bz[i][1] = 0;
			}

	fftw_execute(fftw_plan_dft_3d(n, n, n, Bx, Bx, -1, FFTW_ESTIMATE));
	fftw_execute(fftw_plan_dft_3d(n, n, n, By, By, -1, FFTW_ESTIMATE));
	fftw_execute(fftw_plan_dft_3d(n, n, n, Bz, Bz, -1, FFTW_ESTIMATE));

	// possible discrete wavenumbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	double precision = 1e-12;

	int numNotOrthogonal = 0;
	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				int i = ix * n * n + iy * n + iz;
				// scalar product of B(k) and k has to be zero
				double bdk = Bx[i][0] * K[ix] + By[i][0] * K[iy]
						+ Bz[i][0] * K[iz];
				if (abs(bdk) > precision)
					numNotOrthogonal += 1;
			}

	EXPECT_EQ(0, numNotOrthogonal);

	fftw_free(Bx);
	fftw_free(By);
	fftw_free(Bz);
}

TEST(testTurbulentMagneticField, TurbulentRange) {
	// B(k) = for k < kmin and k > kmax

	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();
	size_t n = B.n;

	// transform back into configuration space
	fftw_complex *Bx, *By, *Bz;
	Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
	By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
	Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);

	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				int i = ix * n * n + iy * n + iz;
				Bx[i][0] = B.field[ix][iy][iz].x();
				Bx[i][1] = 0;
				By[i][0] = B.field[ix][iy][iz].y();
				By[i][1] = 0;
				Bz[i][0] = B.field[ix][iy][iz].z();
				Bz[i][1] = 0;
			}

	fftw_execute(fftw_plan_dft_3d(n, n, n, Bx, Bx, -1, FFTW_ESTIMATE));
	fftw_execute(fftw_plan_dft_3d(n, n, n, By, By, -1, FFTW_ESTIMATE));
	fftw_execute(fftw_plan_dft_3d(n, n, n, Bz, Bz, -1, FFTW_ESTIMATE));

	// possible discrete wavenumbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	double precision = 1e-12;

	int numOutOfRangeNotZero = 0;
	for (size_t ix = 0; ix < n; ix++)
		for (size_t iy = 0; iy < n; iy++)
			for (size_t iz = 0; iz < n; iz++) {
				double k = Vector3(K[ix], K[iy], K[iz]).mag();
				int i = ix * n * n + iy * n + iz;
				Vector3 b(Bx[i][0], Bz[i][0], Bz[i][0]);
				if ((k < 1. / B.lMax) || (k > 1. / B.lMin))
					if (b.mag() > precision)
						numOutOfRangeNotZero += 1;
			}

	EXPECT_EQ(0, numOutOfRangeNotZero);

	fftw_free(Bx);
	fftw_free(By);
	fftw_free(Bz);
}

//TEST(testTurbulentMagneticField, Output) {
//	// Output
//
//	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
//	B.initialize();
//	size_t n = B.n;
//
//	for (size_t i = 0; i < n * n * n; i++)
//		std::cout << B.field[i].x() << ", " << B.field[i].y() << ", "
//				<< B.field[i].z() << std::endl;
//}

TEST(testTurbulentMagneticField, PeriodicBoundaries) {
	// B(x+a*n) = B(x)

	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();
	size_t n = B.n;

	Vector3 position(1.1,2.1,3.1);
	Vector3 b = B.getField(position);
	Vector3 b1 = B.getField(position+Vector3(n,0,0));
	Vector3 b2 = B.getField(position+Vector3(0,2*n,0));
	Vector3 b3 = B.getField(position+Vector3(n,0,3*n));

	EXPECT_FLOAT_EQ(b.x(), b1.x());
	EXPECT_FLOAT_EQ(b.y(), b1.y());
	EXPECT_FLOAT_EQ(b.z(), b1.z());

	EXPECT_FLOAT_EQ(b.x(), b2.x());
	EXPECT_FLOAT_EQ(b.y(), b2.y());
	EXPECT_FLOAT_EQ(b.z(), b2.z());

	EXPECT_FLOAT_EQ(b.x(), b3.x());
	EXPECT_FLOAT_EQ(b.y(), b3.y());
	EXPECT_FLOAT_EQ(b.z(), b3.z());
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

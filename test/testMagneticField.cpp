#include "gtest/gtest.h"
#include "mpc/MagneticField.h"
#include "mpc/TurbulentMagneticField.h"
#include "mpc/Vector3.h"
#include "fftw3.h"
#include <iostream>

namespace mpc {

TEST(testTurbulentMagneticField, Bmean) {
	// turbulent field on a grid with 64^3 nodes with a spacing of 1 and Brms=1
	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();

	Vector3 b(0, 0, 0);
	for (size_t i = 0; i < B.field.size(); i++)
		b += B.field[i];
	b /= B.field.size();

	double precision = 1e-12;
	EXPECT_NEAR(b.x(), 0, precision);
	EXPECT_NEAR(b.y(), 0, precision);
	EXPECT_NEAR(b.z(), 0, precision);
}

TEST(testTurbulentMagneticField, Brms) {
	// turbulent field on a grid with 64^3 nodes with a spacing of 1 and Brms=1
	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
	B.initialize();

	double brms = 0;
	for (size_t i = 0; i < B.field.size(); i++)
		brms += B.field[i].mag2();
	brms = sqrt(brms / B.field.size());

	double precision = 1e-12;
	EXPECT_NEAR(brms, 1, precision);
}

class TurbulentFieldTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		// turbulent field on a grid with 64^3 nodes with a spacing of 1 and Brms=1
		TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
		B.initialize();
		size_t n = B.n;

		// transform back into configuration space
		fftw_complex *Bx, *By, *Bz;
		Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
		By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
		Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);

		for (size_t i = 0; i < B.field.size(); i++) {
			Bx[i][0] = B.field[i].x();
			Bx[i][1] = 0;
			By[i][0] = B.field[i].y();
			By[i][1] = 0;
			Bz[i][0] = B.field[i].z();
			Bz[i][1] = 0;
		}

		fftw_plan p;
		p = fftw_plan_dft_3d(n, n, n, Bx, Bx, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		p = fftw_plan_dft_3d(n, n, n, By, By, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		p = fftw_plan_dft_3d(n, n, n, Bz, Bz, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);

		// possible discrete wavenumbers
		double K[n];
		for (int i = 0; i < n / 2; i++)
			K[i] = (float) i / n;
		for (int i = n / 2; i < n; i++)
			K[i] = (float) i / n - 1;
	}
};

//TEST(testTurbulentMagneticField, turbulentRange) {
//	double precision = 1e-12;
//
//	int numOutOfRangeNotZero = 0;
//	for (size_t ix = 0; ix < n; ix++)
//		for (size_t iy = 0; iy < n; iy++)
//			for (size_t iz = 0; iz < n; iz++) {
//				double k = Vector3(K[ix], K[iy], K[iz]).mag();
//				int i = ix * n * n + iy * n + iz;
//				Vector3 b(Bx[i][0], Bz[i][0], Bz[i][0]);
//
//				if ((k < 1. / B.lMax) || (k > 1. / B.lMin))
//					if (b.mag() > precision)
//						numOutOfRangeNotZero += 1;
//			}
//	EXPECT_EQ(0, numOutOfRangeNotZero);
//}

//TEST(testTurbulentMagneticField, mean) {
//	double precision = 1e-12;
//
//	// turbulent field on a grid with 64^3 nodes with a spacing of 1 and Brms=1
//	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
//	B.initialize();
//	size_t n = B.n;
//
//	// transform back into configuration space
//	fftw_complex *Bx, *By, *Bz;
//	Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
//	By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
//	Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
//
//	for (size_t i = 0; i < B.field.size(); i++) {
//		Bx[i][0] = B.field[i].x();
//		Bx[i][1] = 0;
//		By[i][0] = B.field[i].y();
//		By[i][1] = 0;
//		Bz[i][0] = B.field[i].z();
//		Bz[i][1] = 0;
//	}
//
//	fftw_plan p;
//	p = fftw_plan_dft_3d(n, n, n, Bx, Bx, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	p = fftw_plan_dft_3d(n, n, n, By, By, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	p = fftw_plan_dft_3d(n, n, n, Bz, Bz, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	fftw_destroy_plan(p);
//
//	// possible discrete wavenumbers
//	double K[n];
//	for (int i = 0; i < n / 2; i++)
//		K[i] = (float) i / n;
//	for (int i = n / 2; i < n; i++)
//		K[i] = (float) i / n - 1;
//
//	// <B> = 0 ----------------------------------------------------------------
//	Vector3 b(0, 0, 0);
//	for (size_t i = 0; i < B.field.size(); i++)
//		b += B.field[i];
//	b /= B.field.size();
//	EXPECT_NEAR(b.x(), 0, precision);
//	EXPECT_NEAR(b.y(), 0, precision);
//	EXPECT_NEAR(b.z(), 0, precision);
//
//	// <B^2> = 1 --------------------------------------------------------------
//	double brms = 0;
//	for (size_t i = 0; i < B.field.size(); i++)
//		brms += B.field[i].mag2();
//	brms = sqrt(brms / B.field.size());
//	EXPECT_NEAR(brms, 1, precision);
//
//	// solenoidity ------------------------------------------------------------
//	int numNotOrthogonal = 0;
////	for (size_t ix = 0; ix < n; ix++)
////		for (size_t iy = 0; iy < n; iy++)
////			for (size_t iz = 0; iz < n; iz++) {
////				Vector3 k(K[ix], K[iy], K[iz]);
////				int i = ix * n * n + iy * n + iz;
////				Vector3 b(Bx[i][0], Bz[i][0], Bz[i][0]);
////				if (b.dot(k) > precision)
////					numNotOrthogonal += 1;
////			}
////	EXPECT_EQ(0, numNotOrthogonal);
//
//// turbulent range --------------------------------------------------------
//	int numOutOfRangeNotZero = 0;
//	for (size_t ix = 0; ix < n; ix++)
//		for (size_t iy = 0; iy < n; iy++)
//			for (size_t iz = 0; iz < n; iz++) {
//				double k = Vector3(K[ix], K[iy], K[iz]).mag();
//				int i = ix * n * n + iy * n + iz;
//				Vector3 b(Bx[i][0], Bz[i][0], Bz[i][0]);
//
//				if ((k < 1. / B.lMax) || (k > 1. / B.lMin))
//					if (b.mag() > precision)
//						numOutOfRangeNotZero += 1;
//			}
//	EXPECT_EQ(0, numOutOfRangeNotZero);
//}
//
//TEST(testTurbulentMagneticField, mean) {
//	double precision = 1e-12;
//
//	// turbulent field on a grid with 64^3 nodes with a spacing of 1 and Brms=1
//	TurbulentMagneticField B(Vector3(0, 0, 0), 64, 1, 1, -11. / 3., 2, 8);
//	B.initialize();
//	size_t n = B.n;
//
//	// transform back into configuration space
//	fftw_complex *Bx, *By, *Bz;
//	Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
//	By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
//	Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * n * n);
//
//	for (size_t i = 0; i < B.field.size(); i++) {
//		Bx[i][0] = B.field[i].x();
//		Bx[i][1] = 0;
//		By[i][0] = B.field[i].y();
//		By[i][1] = 0;
//		Bz[i][0] = B.field[i].z();
//		Bz[i][1] = 0;
//	}
//
//	fftw_plan p;
//	p = fftw_plan_dft_3d(n, n, n, Bx, Bx, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	p = fftw_plan_dft_3d(n, n, n, By, By, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	p = fftw_plan_dft_3d(n, n, n, Bz, Bz, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	fftw_destroy_plan(p);
//
//	// possible discrete wavenumbers
//	double K[n];
//	for (int i = 0; i < n / 2; i++)
//		K[i] = (float) i / n;
//	for (int i = n / 2; i < n; i++)
//		K[i] = (float) i / n - 1;
//
//	// <B> = 0 ----------------------------------------------------------------
//	Vector3 b(0, 0, 0);
//	for (size_t i = 0; i < B.field.size(); i++)
//		b += B.field[i];
//	b /= B.field.size();
//	EXPECT_NEAR(b.x(), 0, precision);
//	EXPECT_NEAR(b.y(), 0, precision);
//	EXPECT_NEAR(b.z(), 0, precision);
//
//	// <B^2> = 1 --------------------------------------------------------------
//	double brms = 0;
//	for (size_t i = 0; i < B.field.size(); i++)
//		brms += B.field[i].mag2();
//	brms = sqrt(brms / B.field.size());
//	EXPECT_NEAR(brms, 1, precision);
//
//	// solenoidity ------------------------------------------------------------
//	int numNotOrthogonal = 0;
////	for (size_t ix = 0; ix < n; ix++)
////		for (size_t iy = 0; iy < n; iy++)
////			for (size_t iz = 0; iz < n; iz++) {
////				Vector3 k(K[ix], K[iy], K[iz]);
////				int i = ix * n * n + iy * n + iz;
////				Vector3 b(Bx[i][0], Bz[i][0], Bz[i][0]);
////				if (b.dot(k) > precision)
////					numNotOrthogonal += 1;
////			}
////	EXPECT_EQ(0, numNotOrthogonal);
//
//// turbulent range --------------------------------------------------------
//	int numOutOfRangeNotZero = 0;
//	for (size_t ix = 0; ix < n; ix++)
//		for (size_t iy = 0; iy < n; iy++)
//			for (size_t iz = 0; iz < n; iz++) {
//				double k = Vector3(K[ix], K[iy], K[iz]).mag();
//				int i = ix * n * n + iy * n + iz;
//				Vector3 b(Bx[i][0], Bz[i][0], Bz[i][0]);
//
//				if ((k < 1. / B.lMax) || (k > 1. / B.lMin))
//					if (b.mag() > precision)
//						numOutOfRangeNotZero += 1;
//			}
//	EXPECT_EQ(0, numOutOfRangeNotZero);
//}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/Units.h"
#include "mpc/Common.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(testSPHMagneticField, simpleTest) {
	// Tests if a direct SPH field can be constructed and prints RMS and mean field strength
	Vector3d origin(80 * Mpc);
	double size = 40 * Mpc;

	// gadget::DirectField may throw if queried for magnetic fields on the borders
	Vector3d safeOrigin = origin - Vector3d(1 * kpc);
	double safeSize = size + 2 * kpc;

	SPHMagneticField B(safeOrigin, safeSize, 20, getDataPath("SPH/mhd_z.db").c_str());

	int n = 64;
	double spacing = 40 * Mpc / (n - 1);
	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = B.getField(origin + Vector3d(ix, iy, iz) * spacing);
				brms += b.getMag2();
				bmean += b;
			}

	brms = sqrt(brms / n / n / n);
	bmean /= n * n * n;

	std::cout << "Mean B-Field: " << bmean / nG << " nG" << std::endl;
	std::cout << "RMS B-Field: " << brms / nG << " nG" << std::endl;
}

TEST(testSPHMagneticFieldGrid, simpleTest) {
	// Tests if a sampled SPH field can be constructed and prints RMS and mean field strength
	Vector3d origin(80 * Mpc);
	size_t n = 64;
	double size = 40 * Mpc;
	double spacing = size / (n - 1);

	// gadget::SampledField may throw if queried for magnetic fields on the borders
	Vector3d safeOrigin = origin - Vector3d(1 * kpc);
	double safeSize = size + 2 * kpc;

	SPHMagneticFieldGrid B(safeOrigin, safeSize, n, getDataPath("SPH/mhd_z.db").c_str());

	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = B.getField(origin + Vector3d(ix, iy, iz) * spacing);
				brms += b.getMag2();
				bmean += b;
			}

	brms = sqrt(brms / n / n / n);
	bmean /= n * n * n;

	std::cout << "Mean B-Field: " << bmean / nG << " nG" << std::endl;
	std::cout << "RMS B-Field: " << brms / nG << " nG" << std::endl;
}

TEST(testSPHTurbulentMagneticField, modulatedField) {
	// Test for correct rms and mean strength
	Vector3d origin(80 * Mpc);
	size_t n = 64;

	SPHTurbulentMagneticField B(origin, 40 * Mpc, n);
	double spacing = B.getGridSpacing();
	B.setTurbulenceProperties(2 * spacing, 8 * spacing, -11./3);
	B.initialize();
	B.modulate(getDataPath("SPH/mhd_z.db").c_str(), 2./3);
	B.normalize(1. / B.getRMSFieldStrengthInSphere(Vector3d(120 * Mpc), 105 * Mpc));

	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = B.getField(origin + Vector3d(ix, iy, iz) * spacing);
				brms += b.getMag2();
				bmean += b;
			}

	brms = sqrt(brms / n / n / n);
	bmean /= n * n * n;

	EXPECT_NEAR(brms, 1, 1e-7);
	EXPECT_NEAR(bmean.x, 0, 5e-3); // compatible with 0 within 0.5%
	EXPECT_NEAR(bmean.y, 0, 5e-3);
	EXPECT_NEAR(bmean.z, 0, 5e-3);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/magneticField/SPHTurbulentMagneticFieldGrid.h"
#include "mpc/Units.h"
#include "mpc/Common.h"

#include "gtest/gtest.h"

namespace mpc {

//TEST(testSPHMagneticField, simpleTest) {
//	// Tests if a direct SPH field can be constructed and prints RMS and mean field strength
//	Vector3d origin(80 * Mpc);
//	double size = 40 * Mpc;
//
//	// gadget::DirectField may throw if queried for magnetic fields on the borders
//	Vector3d safeOrigin = origin - Vector3d(1 * kpc);
//	double safeSize = size + 2 * kpc;
//
//	SPHMagneticField B(safeOrigin, safeSize, 20,
//			getDataPath("SPH/mhd_z.db").c_str());
//
//	int n = 64;
//	double spacing = 40 * Mpc / (n - 1);
//	double brms = 0;
//	Vector3d bmean(0, 0, 0);
//	for (int ix = 0; ix < n; ix++)
//		for (int iy = 0; iy < n; iy++)
//			for (int iz = 0; iz < n; iz++) {
//				Vector3d b = B.getField(
//						origin + Vector3d(ix, iy, iz) * spacing);
//				brms += b.getMag2();
//				bmean += b;
//			}
//
//	brms = sqrt(brms / n / n / n);
//	bmean /= n * n * n;
//
//	std::cout << "Mean B-Field: " << bmean / nG << " nG" << std::endl;
//	std::cout << "RMS B-Field: " << brms / nG << " nG" << std::endl;
//}
//
//TEST(testSPHMagneticFieldGrid, simpleTest) {
//	// Tests if a sampled SPH field can be constructed and prints RMS and mean field strength
//	Vector3d origin(80 * Mpc);
//	size_t n = 64;
//	double size = 40 * Mpc;
//	double spacing = size / (n - 1);
//
//	// gadget::SampledField may throw if queried for magnetic fields on the borders
//	Vector3d safeOrigin = origin - Vector3d(1 * kpc);
//	double safeSize = size + 2 * kpc;
//
//	SPHMagneticFieldGrid B(safeOrigin, safeSize, n,
//			getDataPath("SPH/mhd_z.db").c_str());
//
//	double brms = 0;
//	Vector3d bmean(0, 0, 0);
//	for (int ix = 0; ix < n; ix++)
//		for (int iy = 0; iy < n; iy++)
//			for (int iz = 0; iz < n; iz++) {
//				Vector3d b = B.getField(
//						origin + Vector3d(ix, iy, iz) * spacing);
//				brms += b.getMag2();
//				bmean += b;
//			}
//
//	brms = sqrt(brms / n / n / n);
//	bmean /= n * n * n;
//
//	std::cout << "Mean B-Field: " << bmean / nG << " nG" << std::endl;
//	std::cout << "RMS B-Field: " << brms / nG << " nG" << std::endl;
//}

TEST(testSPHTurbulentMagneticField, simpleTest) {
	// Tests if a sampled SPH field can be constructed and prints RMS and mean field strength
	Vector3d origin(100 * Mpc);
	size_t n = 64;
	double size = 40 * Mpc;
	double spacing = size / n;

	SPHTurbulentMagneticFieldGrid B(origin, size, n);
	B.setTurbulenceProperties(2 * spacing, 8 * spacing, -11. / 3);
	B.initialize();
	B.normalize(1 / B.getRMSFieldStrength());

	std::cout << "modulating" << std::endl;
	B.setDensityField(getDataPath("SPH/mhd_z.db").c_str(), Vector3d(99 * Mpc),
			42 * Mpc, 20);
	B.setModulation(getDataPath("SPH/miniati-profile.txt"));

	std::cout << "sampling" << std::endl;
	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = B.getField(
						origin + Vector3d(ix, iy, iz) * spacing);
				brms += b.getMag2();
				bmean += b;
			}

	brms = sqrt(brms / n / n / n);
	bmean /= n * n * n;

	std::cout << "Mean B-Field: " << bmean / nG << " nG" << std::endl;
	std::cout << "RMS B-Field: " << brms / nG << " nG" << std::endl;
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

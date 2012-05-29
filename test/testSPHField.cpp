#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/Units.h"
#include "mpc/Common.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(testSPHMagneticField, simpleTest) {
	// Tests if a direct SPH field can be constructed and prints RMS and mean field strength
	Vector3d origin(80 * Mpc);
	int n = 64;
	double size = 40 * Mpc;
	double spacing = size / n;

	SPHMagneticField bField(origin, 40 * Mpc, 20, getDataPath("SPH/mhd_z.db").c_str());

	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = bField.getField(origin + Vector3d(ix, iy, iz) * spacing);
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
	Vector3d margin(1 * kpc);
	int n = 64;
	double size = 40 * Mpc;
	double spacing = size / n;

	// gadget::SampledField throws errors if queried for magnetic fields on the borders
	Vector3d safeOrigin(80.001 * Mpc);
	double safeSpacing = (size - 2 * kpc) / n;

	SPHMagneticFieldGrid bField(origin, 40 * Mpc, n, getDataPath("SPH/mhd_z.db").c_str());

	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = bField.getField(safeOrigin + Vector3d(ix, iy, iz) * safeSpacing);
				brms += b.getMag2();
				bmean += b;
			}

	brms = sqrt(brms / n / n / n);
	bmean /= n * n * n;

	std::cout << "Mean B-Field: " << bmean / nG << " nG" << std::endl;
	std::cout << "RMS B-Field: " << brms / nG << " nG" << std::endl;
}

TEST(testSPHTurbulentMagneticField, modulatedField) {
	// Test SPH modulated turbulent field for correct RMS strength: <B^2> = Brms^2 and mean
	Vector3d origin(80 * Mpc);
	double size = 40 * Mpc; // subvolume of the SPH field
	int n = 64;
	double spacing = size / n;
	double lMin = 2 * spacing;
	double lMax = 16 * spacing;
	double Brms = 1;
	double alpha = -11./3;

	SPHTurbulentMagneticField bField(origin, n, spacing);
	bField.initialize(lMin, lMax, Brms, alpha);
	bField.modulate(getDataPath("SPH/mhd_z.db").c_str());

	double brms = 0;
	Vector3d bmean(0, 0, 0);
	for (int ix = 0; ix < n; ix++)
		for (int iy = 0; iy < n; iy++)
			for (int iz = 0; iz < n; iz++) {
				Vector3d b = bField.getField(origin + Vector3d(ix, iy, iz) * spacing);
				brms += b.getMag2();
				bmean += b;
			}

	brms = sqrt(brms / n / n / n);
	bmean /= n * n * n;

	EXPECT_NEAR(brms, 1, 1e-7);
	EXPECT_NEAR(bmean.x, 0, 1e-2); // compatible with 0 within 1%
	EXPECT_NEAR(bmean.y, 0, 1e-2);
	EXPECT_NEAR(bmean.z, 0, 1e-2);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

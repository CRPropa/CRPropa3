#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

namespace mpc {

void SPHTurbulentMagneticField::modulate(const std::string filename) {
	// create SPH Field to obtain baryon density
	std::cout << "mpc::SPHTurbulentMagneticField: Loading SPH-field." << std::endl;
	Vector3d safeOrigin = origin - Vector3d(1,1,1) * kpc;
	double safeSize = spacing * samples + 2 * kpc;
	SPHMagneticField sph(safeOrigin, safeSize, 32, filename);

	// modulate and calculate renormalization
	std::cout << "mpc::SPHTurbulentMagneticField: Modulate turbulent field." << std::endl;
	double sumB2 = 0;
	int n = 0;
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++) {
				Vector3d pos = Vector3d(ix, iy, iz) * spacing + origin;
				double rho = sph.getRho(pos);
				Vector3f &b = get(ix, iy, iz);
				b *= pow(rho, 2. / 3);
				double dist = (pos / Mpc - Vector3d(120, 120, 120)).getMag();
				if (dist < 110) {
					sumB2 += b.getMag2();
					n++;
				}
			}

	// renormalize
	double norm = Brms / sqrt(sumB2 / n);
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++) {
				Vector3f &b = get(ix, iy, iz);
				b *= norm;
			}
}

} // namespace mpc

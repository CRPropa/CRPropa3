#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

namespace mpc {

void SPHTurbulentMagneticField::modulate(std::string filename, double weight) {
	// SPHMagneticField fails on borders, choose minimally larger size
	Vector3d safeOrigin = origin - Vector3d(kpc);
	double safeSize = size + 2 * kpc;
	SPHMagneticField sph(safeOrigin, safeSize, 20, filename);

	// modulate
	for (int ix = 0; ix < samples; ix++) {
		for (int iy = 0; iy < samples; iy++) {
			for (int iz = 0; iz < samples; iz++) {
				Vector3d pos = origin + Vector3d(ix, iy, iz) * spacing;
				double rho = sph.getRho(pos);
				get(ix, iy, iz) *= pow(rho, weight);
			}
		}
	}
}

} // namespace mpc

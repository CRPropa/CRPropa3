#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/magneticField/SPHMagneticField.h"

namespace mpc {

void SPHTurbulentMagneticField::modulate(const std::string filename) {
	// create SPH Field to obtain baryon density
	SPHMagneticField sphField(origin, spacing * samples, samples, filename);

	// modulate and calculate renormalization
	double norm = 0;
	for (int ix; ix < samples; ix++)
		for (int iy; iy < samples; iy++)
			for (int iz; iz < samples; iz++) {
				double rho = sphField.getRho(
						Vector3d(ix, iy, iz) * spacing + origin);
				Vector3f &b = get(ix, iy, iz);
				b *= rho;
				norm += b.mag2();
			}

	// renormalize
	norm = Brms / sqrt(norm / pow(samples, 3));
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++) {
				Vector3f &b = get(ix, iy, iz);
				b *= norm;
			}
}

} // namespace mpc

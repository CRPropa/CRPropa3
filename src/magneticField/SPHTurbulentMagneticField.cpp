#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

namespace mpc {

void SPHTurbulentMagneticField::modulate(std::string filename) {
	// create SPH Field to obtain baryon density
	std::cout << "mpc::SPHTurbulentMagneticField: Loading SPH-field." << std::endl;
	SPHMagneticField sph(origin, spacing * samples, 20, filename);

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
				double dist = (pos - Vector3d(120 * Mpc)).getMag();
				if (dist < 110 * Mpc) {
					sumB2 += b.getMag2();
					n++;
				}
			}

	// renormalize
	std::cout << "mpc::SPHTurbulentMagneticField: Normalize turbulent field." << std::endl;
	double norm = Brms / sqrt(sumB2 / n);
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++) {
				Vector3f &b = get(ix, iy, iz);
				b *= norm;
			}
}

} // namespace mpc

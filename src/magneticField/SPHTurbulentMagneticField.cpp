#include "mpc/magneticField/SPHTurbulentMagneticField.h"
#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

namespace mpc {

void SPHTurbulentMagneticField::modulate(std::string filename, double weight = 2./3) {
	// create SPH Field to obtain baryon density
	std::cout << "mpc::SPHTurbulentMagneticField: Loading SPH-field: "
			<< filename << std::endl;
	// SPHMagneticField may fail on borders, choose minimally larger size
	Vector3d safeOrigin = origin - Vector3d(1.);
	double safeSize = size + 2;
	SPHMagneticField sph(safeOrigin, safeSize, 20, filename);

	// modulate and calculate renormalization
	std::cout << "mpc::SPHTurbulentMagneticField: Modulate turbulent field."
			<< std::endl;
	double sumB2 = 0;
	int numB2 = 0;
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++) {
				Vector3d pos = origin + Vector3d(ix, iy, iz) * spacing;
				double rho = sph.getRho(pos);
				Vector3f &b = get(ix, iy, iz);
				b *= pow(rho, weight);
				double distance2center = (pos - Vector3d(120 * Mpc)).getMag();
				if (distance2center < 110 * Mpc) {
					sumB2 += b.getMag2();
					numB2++;
				}
			}

	// renormalize
	std::cout << "mpc::SPHTurbulentMagneticField: Normalize turbulent field."
			<< std::endl;
	double norm = Brms / sqrt(sumB2 / numB2);
	Brms = 0;
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++) {
				Vector3f &b = get(ix, iy, iz);
				b *= norm;
				Brms += b.getMag2();
			}
}

} // namespace mpc

#include "mpc/magneticField/SPHTurbulentMagneticFieldGrid.h"
#include "mpc/Common.h"
#include "mpc/Units.h"

namespace mpc {

void SPHTurbulentMagneticFieldGrid::setDensityField(std::string filename,
		Vector3d origin, double size, size_t bins) {
	sphField.reset(new SPHMagneticField(origin, size, bins, filename));
}

void SPHTurbulentMagneticFieldGrid::setModulation(std::string filename) {
	// load modulation profile
	std::ifstream infile(filename.c_str());
	if (!infile.good())
		throw std::runtime_error("mpc: could not open file " + filename);
	for (int i = 0; i < 67; i++) {
		infile >> tabRho[i] >> tabB[i];
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

Vector3d SPHTurbulentMagneticFieldGrid::getField(
		const Vector3d &position) const {
	double rho = sphField->getRho(position); // density in [kg/m^3]
	if (rho <= std::numeric_limits<double>::min())
		return Vector3d(0.);
	double x = log10(rho / 1000); // log10(density in [g/cm^3])
	double m; // strength in [G]
	if (x <= -32)
		m = 206.404 + 12.9872 * x + 0.193117 * x * x;
	else if (x >= -28)
		m = -37.3625 - 2.74343 * x - 0.0585562 * x * x;
	else
		m = interpolate(x, tabRho, tabB);
	return MagneticFieldGrid::getField(position) * pow(10, m) / 10000; // B-field in [T]
}

double SPHTurbulentMagneticFieldGrid::getRho(const Vector3d &position) const {
	return sphField->getRho(position);
}

} // namespace mpc

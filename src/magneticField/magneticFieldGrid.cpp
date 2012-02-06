#include "mpc/magneticField/magneticFieldGrid.hpp"

namespace mpc {

void periodicClamp(double pos, int n, int& lower, int& upper, double& fracLower,
		double& fracUpper) {
	// closest lower and upper neighbour in a periodically continued grid
	lower = ((int(pos) % n) + n) % n;
	upper = (lower + 1) % n;
	fracLower = pos - floor(pos);
	fracUpper = 1 - fracLower;
}

MagneticFieldGrid::MagneticFieldGrid(size_t n, double spacing,
		Vector3 origin) {
	this->origin = origin;
	this->spacing = spacing;
	this->n = n;
	grid.resize(n);
	for (int ix = 0; ix < n; ix++) {
		grid[ix].resize(n);
		for (int iy = 0; iy < n; iy++) {
			grid[ix][iy].resize(n);
		}
	}
}

Vector3 MagneticFieldGrid::getField(const Vector3 &position) const {
	Vector3 r = (position - origin) / spacing;
	int ix, iX, iy, iY, iz, iZ;
	double fx, fX, fy, fY, fz, fZ;
	periodicClamp(r.x(), n, ix, iX, fx, fX);
	periodicClamp(r.y(), n, iy, iY, fy, fY);
	periodicClamp(r.z(), n, iz, iZ, fz, fZ);

	// trilinear interpolation
	// check: http://paulbourke.net/miscellaneous/interpolation/
	Vector3 b(0, 0, 0);
	//V000 (1 - x) (1 - y) (1 - z) +
	b += grid[ix][iy][iz] * fX * fY * fZ;
	//V100 x (1 - y) (1 - z) +
	b += grid[iX][iy][iz] * fx * fY * fZ;
	//V010 (1 - x) y (1 - z) +
	b += grid[ix][iY][iz] * fX * fy * fZ;
	//V001 (1 - x) (1 - y) z +
	b += grid[ix][iy][iZ] * fX * fY * fz;
	//V101 x (1 - y) z +
	b += grid[iX][iy][iZ] * fx * fY * fz;
	//V011 (1 - x) y z +
	b += grid[ix][iY][iZ] * fX * fy * fz;
	//V110 x y (1 - z) +
	b += grid[iX][iY][iz] * fx * fy * fZ;
	//V111 x y z
	b += grid[iX][iY][iZ] * fx * fy * fz;

	return b;
}

} // namespace mpc

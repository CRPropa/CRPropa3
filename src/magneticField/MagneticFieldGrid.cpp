#include "mpc/magneticField/MagneticFieldGrid.h"

namespace mpc {

void periodicClamp(double x, int n, int &low, int &high, double &fLow,
		double &fHigh) {
	// closest lower and upper neighbour in a periodically continued grid
	low = ((int(x) % n) + n) % n;
	high = (low + 1) % n;
	fLow = x - floor(x);
	fHigh = 1 - fLow;
}

MagneticFieldGrid::MagneticFieldGrid(Vector3d origin, size_t samples,
		double spacing) :
		origin(origin), samples(samples), spacing(spacing) {
	grid.resize(samples * samples * samples);
}

void MagneticFieldGrid::updateSimulationVolume(const Vector3d &origin,
		double size) {
	this->origin = origin;
	this->spacing = size / (samples - 1);
}

Vector3d MagneticFieldGrid::getGridOrigin() const {
	return origin;
}

size_t MagneticFieldGrid::getGridSamples() const {
	return samples;
}

double MagneticFieldGrid::getGridSpacing() const {
	return spacing;
}

double MagneticFieldGrid::getGridSize() const {
	return samples * spacing;
}

Vector3f &MagneticFieldGrid::get(size_t ix, size_t iy, size_t iz) {
	int i = ix * samples * samples + iy * samples + iz;
	return grid[i];
}

const Vector3f &MagneticFieldGrid::get(size_t ix, size_t iy, size_t iz) const {
	int i = ix * samples * samples + iy * samples + iz;
	return grid[i];
}

Vector3d MagneticFieldGrid::getField(const Vector3d &position) const {
	Vector3d r = (position - origin) / spacing;
	int ix, iX, iy, iY, iz, iZ;
	double fx, fX, fy, fY, fz, fZ;
	periodicClamp(r.x, samples, ix, iX, fx, fX);
	periodicClamp(r.y, samples, iy, iY, fy, fY);
	periodicClamp(r.z, samples, iz, iZ, fz, fZ);

	// trilinear interpolation
	// check: http://paulbourke.net/miscellaneous/interpolation/
	Vector3d b(0, 0, 0);
	//V000 (1 - x) (1 - y) (1 - z) +
	b += get(ix, iy, iz) * fX * fY * fZ;
	//V100 x (1 - y) (1 - z) +
	b += get(iX, iy, iz) * fx * fY * fZ;
	//V010 (1 - x) y (1 - z) +
	b += get(ix, iY, iz) * fX * fy * fZ;
	//V001 (1 - x) (1 - y) z +
	b += get(ix, iy, iZ) * fX * fY * fz;
	//V101 x (1 - y) z +
	b += get(iX, iy, iZ) * fx * fY * fz;
	//V011 (1 - x) y z +
	b += get(ix, iY, iZ) * fX * fy * fz;
	//V110 x y (1 - z) +
	b += get(iX, iY, iz) * fx * fy * fZ;
	//V111 x y z
	b += get(iX, iY, iZ) * fx * fy * fz;

	return b;
}

void MagneticFieldGrid::setGridOrigin(const Vector3d& origin) {
	this->origin = origin;
}

void MagneticFieldGrid::setGridSpacing(const double spacing) {
	this->spacing = spacing;
}

} // namespace mpc

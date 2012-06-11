#include "mpc/magneticField/MagneticFieldGrid.h"

#include <stdexcept>
#include <fstream>

namespace mpc {

void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ( ( int(floor(x)) % n ) + n ) % n;
	hi = (lo + 1) % n;
}

MagneticFieldGrid::MagneticFieldGrid(Vector3d origin, double size,
		size_t samples) {
	this->origin = origin;
	this->size = size;
	this->samples = samples;
	this->spacing = size / samples;
	grid.resize(samples * samples * samples);
}

void MagneticFieldGrid::updateSimulationVolume(const Vector3d &origin,
		double size) {
	this->origin = origin;
	this->size = size;
	this->spacing = size / samples;
}

void MagneticFieldGrid::normalize(double norm) {
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++)
				get(ix, iy, iz) *= norm;
}

void MagneticFieldGrid::modulateWithDensityField(std::string filename,
		double exp) {
	std::ifstream infile(filename.c_str(), std::ios::binary);
	if (!infile)
		throw std::runtime_error("mpc::MagneticFieldGrid: File not found");
	double sumB2 = 0;
	int numB2 = 0;
	float rho;
	for (int ix = 0; ix < samples + 1; ix++) {
		for (int iy = 0; iy < samples + 1; iy++) {
			for (int iz = 0; iz < samples + 1; iz++) {
				infile.read((char*) &rho, sizeof(float));
				if ((ix == samples) or (iy == samples) or (iz == samples))
					continue; // skip last grid points in each direction
				Vector3d pos = origin + Vector3d(ix, iy, iz) * spacing;
				Vector3f &b = get(ix, iy, iz);
				b *= pow(rho, exp);
			}
		}
	}
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
	return size;
}

double MagneticFieldGrid::getRMSFieldStrength() const {
	double sumB2 = 0;
	for (int ix = 0; ix < samples; ix++)
		for (int iy = 0; iy < samples; iy++)
			for (int iz = 0; iz < samples; iz++)
				sumB2 += get(ix, iy, iz).getMag2();
	return sqrt(sumB2 / samples / samples / samples);
}

double MagneticFieldGrid::getRMSFieldStrengthInSphere(Vector3d center,
		double radius) const {
	int numB2 = 0;
	double sumB2 = 0;
	for (int ix = 0; ix < samples; ix++) {
		for (int iy = 0; iy < samples; iy++) {
			for (int iz = 0; iz < samples; iz++) {
				Vector3d position = origin + Vector3d(ix, iy, iz) * spacing;
				if (position.getDistanceTo(center) < radius) {
					sumB2 += get(ix, iy, iz).getMag2();
					numB2++;
				}
			}
		}
	}
	return sqrt(sumB2 / numB2);
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
	// position on a unit grid
	Vector3d r = (position - origin) / spacing;

	// indices of lower and upper neighbors
	int ix, iX, iy, iY, iz, iZ;
	periodicClamp(r.x, samples, ix, iX);
	periodicClamp(r.y, samples, iy, iY);
	periodicClamp(r.z, samples, iz, iZ);

	// linear fraction to lower and upper neighbor
	double fx = r.x - floor(r.x);
	double fX = 1 - fx;
	double fy = r.y - floor(r.y);
	double fY = 1 - fy;
	double fz = r.z - floor(r.z);
	double fZ = 1 - fz;

	// trilinear interpolation
	// check: http://paulbourke.net/miscellaneous/interpolation/
	Vector3d b(0.);
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

} // namespace mpc

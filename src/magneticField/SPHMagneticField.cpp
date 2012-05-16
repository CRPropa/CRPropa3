#include "mpc/magneticField/SPHMagneticField.h"
#include "mpc/Units.h"

#include "kiss/logger.h"

namespace mpc {

SPHMagneticField::SPHMagneticField(Vector3d origin, double size,
		size_t gridSize, const std::string filename) :
		field(gridSize) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z) / kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticField::SPHMagneticField(size_t gridSize, const std::string filename) :
		field(gridSize) {
	database.open(filename);
}

Vector3d SPHMagneticField::getField(const Vector3d &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y, position.z);
	gadget::Vector3f b = field.getField(r / kpc);
	Vector3d bField = Vector3d(b.x, b.y, b.z) * gauss;
	return bField;
}

double SPHMagneticField::getRho(const Vector3d& position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y, position.z);
	size_t overlaps = 0;
	return (double) field.getRho(r / kpc, overlaps);
}

void SPHMagneticField::updateSimulationVolume(const Vector3d &origin,
		double size) {
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z) / kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(Vector3d origin, double size,
		size_t samples, const std::string filename) :
		samples(samples), field(samples), cacheEnabled(false) {
	database.open(filename);
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z) / kpc;
	field.init(v, size / kpc, database);
}

SPHMagneticFieldGrid::SPHMagneticFieldGrid(size_t samples,
		const std::string filename) :
		samples(samples), field(samples) {
	database.open(filename);
}

Vector3d SPHMagneticFieldGrid::getField(const Vector3d &position) const {
	gadget::Vector3f r = gadget::Vector3f(position.x, position.y, position.z);
	gadget::Vector3f b = field.getField(r / kpc);
	Vector3d bField = Vector3d(b.x, b.y, b.z) * gauss;
	return bField;
}

void SPHMagneticFieldGrid::updateSimulationVolume(const Vector3d &origin,
		double size) {
	gadget::Vector3f v = gadget::Vector3f(origin.x, origin.y, origin.z) / kpc;
	float s = size / kpc;
	if (cacheEnabled) {
		std::stringstream path;
		path << cachePrefix << samples << "_" << s << "_" << v.x << "_" << v.y
				<< "_" << v.z << ".cache";
		std::string filename = path.str();
		std::ifstream infile(filename.c_str());
		if (infile) {
			std::cout << "load cache file " << filename << std::endl;
			field.init(v, s);
			field.restore(filename);
		} else {
			field.init(v, s, database);
			path << "." << time(NULL) << clock();
			std::string filename_tmp = path.str();
			field.dump(filename_tmp);
			int ren = rename(filename_tmp.c_str(), filename.c_str());
			if (ren != 0)
				remove(filename_tmp.c_str());
		}
	} else {
		field.init(v, s, database);
	}
}

void SPHMagneticFieldGrid::setCachePrefix(const std::string &prefix) {
	cachePrefix = prefix;
}

void SPHMagneticFieldGrid::setCacheEnabled(bool enabled) {
	cacheEnabled = enabled;
}

} // namespace mpc


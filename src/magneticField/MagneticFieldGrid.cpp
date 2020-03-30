#include "crpropa/magneticField/MagneticFieldGrid.h"

namespace crpropa {

MagneticFieldGrid::MagneticFieldGrid(ref_ptr<VectorGridf> grid) {
	setGrid(grid);
}

void MagneticFieldGrid::setGrid(ref_ptr<VectorGridf> grid) {
	this->grid = grid;
}

ref_ptr<VectorGridf> MagneticFieldGrid::getGrid() {
	return grid;
}

Vector3d MagneticFieldGrid::getField(const Vector3d &pos) const {
	return grid->interpolate(pos);
}

ModulatedMagneticFieldGrid::ModulatedMagneticFieldGrid(ref_ptr<VectorGridf> grid,
		ref_ptr<ScalarGridf> modGrid) {
	grid->setReflective(false);
	modGrid->setReflective(true);
	setGrid(grid);
	setModulationGrid(modGrid);
}

void ModulatedMagneticFieldGrid::setGrid(ref_ptr<VectorGridf> g) {
	grid = g;
}

ref_ptr<VectorGridf> ModulatedMagneticFieldGrid::getGrid() {
	return grid;
}

void ModulatedMagneticFieldGrid::setModulationGrid(ref_ptr<ScalarGridf> g) {
	modGrid = g;
}

ref_ptr<ScalarGridf> ModulatedMagneticFieldGrid::getModulationGrid() {
	return modGrid;
}

void ModulatedMagneticFieldGrid::setReflective(bool gridReflective,
		bool modGridReflective) {
	grid->setReflective(gridReflective);
	modGrid->setReflective(modGridReflective);
}

Vector3d ModulatedMagneticFieldGrid::getField(const Vector3d &pos) const {
	float m = modGrid->interpolate(pos);
	Vector3d b = grid->interpolate(pos);
	return b * m;
}

} // namespace crpropa

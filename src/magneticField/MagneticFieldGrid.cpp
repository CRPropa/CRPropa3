#include "crpropa/magneticField/MagneticFieldGrid.h"

namespace crpropa {

MagneticFieldGrid::MagneticFieldGrid(ref_ptr<Grid3f> grid) {
	setGrid(grid);
}

void MagneticFieldGrid::setGrid(ref_ptr<Grid3f> grid) {
	this->grid = grid;
}

ref_ptr<Grid3f> MagneticFieldGrid::getGrid() {
	return grid;
}

Vector3d MagneticFieldGrid::getField(const Vector3d &pos) const {
	return grid->interpolate(pos);
}

ModulatedMagneticFieldGrid::ModulatedMagneticFieldGrid(ref_ptr<Grid3f> grid,
		ref_ptr<Grid1f> modGrid) {
	grid->setReflective(false);
	modGrid->setReflective(true);
	setGrid(grid);
	setModulationGrid(modGrid);
}

void ModulatedMagneticFieldGrid::setGrid(ref_ptr<Grid3f> g) {
	grid = g;
}

ref_ptr<Grid3f> ModulatedMagneticFieldGrid::getGrid() {
	return grid;
}

void ModulatedMagneticFieldGrid::setModulationGrid(ref_ptr<Grid1f> g) {
	modGrid = g;
}

ref_ptr<Grid1f> ModulatedMagneticFieldGrid::getModulationGrid() {
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

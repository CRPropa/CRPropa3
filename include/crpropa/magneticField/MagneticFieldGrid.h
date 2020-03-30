#ifndef CRPROPA_MAGNETICFIELDGRID_H
#define CRPROPA_MAGNETICFIELDGRID_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class MagneticFieldGrid
 @brief Magnetic field on a periodic (or reflective), cartesian grid with trilinear interpolation.

 This class wraps a VectorGridf to serve as a MagneticField.
 */
class MagneticFieldGrid: public MagneticField {
	ref_ptr<VectorGridf> grid;
public:
	MagneticFieldGrid(ref_ptr<VectorGridf> grid);
	void setGrid(ref_ptr<VectorGridf> grid);
	ref_ptr<VectorGridf> getGrid();
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class ModulatedMagneticFieldGrid
 @brief Modulated magnetic field on a periodic grid.

 This class wraps a VectorGridf to serve as a MagneticField.
 The field is modulated on-the-fly with a ScalarGridf.
 The VectorGridf and ScalarGridf do not need to share the same origin, spacing or size.
 */
class ModulatedMagneticFieldGrid: public MagneticField {
	ref_ptr<VectorGridf> grid;
	ref_ptr<ScalarGridf> modGrid;
public:
	ModulatedMagneticFieldGrid() {
	}
	ModulatedMagneticFieldGrid(ref_ptr<VectorGridf> grid, ref_ptr<ScalarGridf> modGrid);
	void setGrid(ref_ptr<VectorGridf> grid);
	void setModulationGrid(ref_ptr<ScalarGridf> modGrid);
	ref_ptr<VectorGridf> getGrid();
	ref_ptr<ScalarGridf> getModulationGrid();
	void setReflective(bool gridReflective, bool modGridReflective);
	Vector3d getField(const Vector3d &position) const;
};
/** @} */
} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELDGRID_H

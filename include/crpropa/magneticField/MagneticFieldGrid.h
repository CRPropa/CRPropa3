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

 This class wraps a Grid3f to serve as a MagneticField.
 */
class MagneticFieldGrid: public MagneticField {
	ref_ptr<Grid3f> grid;
public:
	/**
	 *Constructor
	 @param grid Grid3f storing the magnetic field vectors
	*/
	MagneticFieldGrid(ref_ptr<Grid3f> grid);
	void setGrid(ref_ptr<Grid3f> grid);
	ref_ptr<Grid3f> getGrid();
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class ModulatedMagneticFieldGrid
 @brief Modulated magnetic field on a periodic grid.

 This class wraps a Grid3f to serve as a MagneticField.
 The field is modulated on-the-fly with a Grid1f.
 The Grid3f and Grid1f do not need to share the same origin, spacing or size.
 */
class ModulatedMagneticFieldGrid: public MagneticField {
	ref_ptr<Grid3f> grid;
	ref_ptr<Grid1f> modGrid;
public:
	ModulatedMagneticFieldGrid() {
	}
	/**
	 *Constructor
	 @param grid 	Grid3f storing the magnetic field vectors
	 @param modGrid Grid1f used to scale the magnetic field strength
	 				B^new_i = B^old_i * scale 
	*/
	ModulatedMagneticFieldGrid(ref_ptr<Grid3f> grid, ref_ptr<Grid1f> modGrid);
	void setGrid(ref_ptr<Grid3f> grid);
	void setModulationGrid(ref_ptr<Grid1f> modGrid);
	ref_ptr<Grid3f> getGrid();
	ref_ptr<Grid1f> getModulationGrid();
	void setReflective(bool gridReflective, bool modGridReflective);
	Vector3d getField(const Vector3d &position) const;
};
/** @} */
} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELDGRID_H

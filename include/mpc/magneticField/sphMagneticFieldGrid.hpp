#ifndef SPHMAGNETICFIELDGRID_HPP_
#define SPHMAGNETICFIELDGRID_HPP_

#include "mpc/magneticField/magneticFieldGrid.hpp"

namespace mpc {

class SPH_MagneticFieldGrid: public MagneticFieldGrid {
public:
	SPH_MagneticFieldGrid(size_t n, double spacing, Vector3 origin,
			std::string filename);
protected:
	void initialize();
	std::string filename;
};

} // namespace mpc

#endif /* SPHMAGNETICFIELDGRID_HPP_ */

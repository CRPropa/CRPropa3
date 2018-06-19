#ifndef CRPROPA_RESTRICTTOREGION_H
#define CRPROPA_RESTRICTTOREGION_H

#include "crpropa/Referenced.h"
#include "crpropa/Candidate.h"
#include "crpropa/Module.h"
#include "crpropa/Geometry.h"

namespace crpropa {
/**
 * \addtogroup Condition
 * @{
 */

/**
 @class RestrictToRegion
 @brief Limit Module to region in simulation.

 */
class RestrictToRegion: public Module {
private:
	ref_ptr<Surface> surface;
	ref_ptr<Module> module;
public:
	RestrictToRegion(Module* _module, Surface* _surface);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_RESTRICTTOREGION_H

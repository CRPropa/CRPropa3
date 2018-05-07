#include "crpropa/Referenced.h"
#include "crpropa/Module.h"
#include "crpropa/Surface.h"

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
	ref_ptr<module> module;
public:
	RestrictToRegion(ref_ptr<Module> _module, ref_ptr<Surface> _surface);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/** @}*/

} // namespace crpropa


#include "crpropa/module/RestrictToRegion.h" 

#include <sstream>

namespace crpropa {

RestrictToRegion::RestrictToRegion(ref_ptr<Module> _module, ref_ptr<Surface> _surface) : module(_module), surface(_surface) { 
};

void RestrictToRegion::process(ref_ptr<Candidate> candidate) const {
	if (surface->distance(candidate->current.getPosition()) <= 0) {
		module->process(candidate);
	}
};

std::string RestrictToRegion::getDescription() const {
	std::stringstream s;
	s << "RestrictToArea:\n"
		<< "  Module: " << module->getDescription() << std::endl
		<< "  Region: " << surface->getDescription() << std::endl;
	return s.str();
};


} // namespace crpropa

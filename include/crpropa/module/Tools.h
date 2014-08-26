#ifndef CRPROPA_MODULETOOLS_H
#define CRPROPA_MODULETOOLS_H

#include "crpropa/Module.h"

#include <vector>

namespace crpropa {

class PerformanceModule: public Module {
private:
	struct _module_info {
		double time;
		ref_ptr<Module> module;
	};

	mutable std::vector<_module_info> modules;
	mutable size_t calls;

public:
	~PerformanceModule();
	void add(Module* module);
	void process(Candidate* candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_MODULETOOLS_H

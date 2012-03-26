#ifndef MPC_MODULE_TOOLS_H_
#define MPC_MODULE_TOOLS_H_

#include "mpc/Module.h"

namespace mpc {

class PerformanceModule: public Module {
private:
	struct _module_info {
		double time;
		ref_ptr<Module> module;
	};

	mutable std::vector<_module_info> modules;
	mutable size_t calls;

public:
	PerformanceModule();
	~PerformanceModule();
	void add(Module *module);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_MODULE_TOOLS_H_ */

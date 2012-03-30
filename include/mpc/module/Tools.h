#ifndef MPC_MODULE_TOOLS_H_
#define MPC_MODULE_TOOLS_H_

#include "mpc/Module.h"
#include "mpc/AssocVector.h"

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

class PropertyStatistics: public Module {
private:
	mutable Loki::AssocVector<std::string, size_t> properties;
	std::string key;

public:
	PropertyStatistics(const std::string &key);
	~PropertyStatistics();
	void process(Candidate *candidate) const;
};

} // namespace mpc

#endif /* MPC_MODULE_TOOLS_H_ */

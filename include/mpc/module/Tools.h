#ifndef MPC_MODULE_TOOLS_H_
#define MPC_MODULE_TOOLS_H_

#include "mpc/Module.h"
#include "mpc/AssocVector.h"

#include <vector>

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
	~PerformanceModule();
	void add(Module* module);
	void process(Candidate* candidate) const;
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

/**
 @class ParticleSelector
 @brief Wraps a module and only executes it for certain particle types
 */
class ParticleSelector: public Module {
private:
	std::vector<int> ids;
	ref_ptr<Module> module;

public:
	ParticleSelector(Module* module);
	ParticleSelector(Module* module, std::vector<int> particleIDs);
	void add(int particleId);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace mpc

#endif /* MPC_MODULE_TOOLS_H_ */

#ifndef CRPROPA_MODULETOOLS_H
#define CRPROPA_MODULETOOLS_H

#include "crpropa/Module.h"
#include "crpropa/EmissionMap.h"

#include <vector>
#include <set>

namespace crpropa {
/**
 * \addtogroup Tools
 * @{
 */

/**
 @class PerformanceModule
 @brief Module to monitor the simulation performance

 Add modules under investigation to this module instead of the ModuleList.
 */
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

/**
  @class ParticleFilter
  @brief Reject Particles not listed in filter.
*/
class ParticleFilter: public AbstractCondition {
	std::set<int> ids;

public:
	ParticleFilter();
	ParticleFilter(const std::set<int> &ids);
	void addId(int id);
	void removeId(int remove);
	std::set<int> &getIds();

	void process(Candidate* candidate) const;
	std::string getDescription() const;
};


/**
  @class EmissionMapFiller
  @brief Fill EmissionMap with source particle state
*/
class EmissionMapFiller: public Module {
	ref_ptr<EmissionMap> emissionMap;
public:
	EmissionMapFiller(EmissionMap *emissionMap);
	void setEmissionMap(EmissionMap *emissionMap);
	void process(Candidate* candidate) const;
	std::string getDescription() const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_MODULETOOLS_H

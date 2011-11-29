#ifndef MODULE_CHAIN_H_
#define MODULE_CHAIN_H_

#include <list>
#include <typeinfo>
#include <assert.h>

#include "mpc/Candidate.h"
#include "mpc/Module.h"

namespace mpc {

struct Priority {
	enum Enum {
		Start = 0,
		BeforeIntegration = 20,
		Integration = 25,
		AfterIntegration = 30,
		BeforeInteraction = 45,
		Interaction = 50,
		AfterInteraction = 55,
		BeforeCommit = 70,
		Commit = 75,
		AfterCommit = 80,
		End = 100
	};
};

class ModuleChain {
public:
	typedef std::pair<size_t, Module *> list_entry_t;
	typedef std::list<list_entry_t> list_t;

	const list_t &getStartModules() const;
	const list_t &getMainModules() const;
	const list_t &getEndModules() const;

	void add(size_t priority, Module *module);
	void apply(Candidate &candidate);
	void clear();

private:
	void check();

	list_t startModules;
	list_t mainModules;
	list_t endModules;
};

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleChain &chain);

#endif /* MODULE_CHAIN_H_ */

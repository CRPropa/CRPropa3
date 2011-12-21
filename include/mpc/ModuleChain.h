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
		BeforePropagation = 20,
		Propagation = 25,
		AfterPropagation = 30,
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

	void add(Module *module, size_t priority);
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	void process(std::vector<Candidate *> &candidates);
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

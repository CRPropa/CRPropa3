#ifndef MPC_MODULE_CHAIN_H_
#define MPC_MODULE_CHAIN_H_

#include <list>
#include <typeinfo>
#include <assert.h>

#include "mpc/Candidate.h"
#include "mpc/Module.h"

namespace mpc {

/**
 @class ModuleChain
 @brief The simulation chain.
 */
class ModuleChain {
public:
	struct Priority {
		enum Enum {
			Start = 0, Propagation = 1, End = 100
		};
	};

	typedef std::pair<size_t, ref_ptr<Module> > list_entry_t;
	typedef std::list<list_entry_t> list_t;

	const list_t &getStartModules() const;
	const list_t &getMainModules() const;
	const list_t &getEndModules() const;

private:
	list_t startModules;
	list_t mainModules;
	list_t endModules;

public:
	void add(size_t priority, Module *module);
	void process(Candidate *candidate);
	void process(std::vector<ref_ptr<Candidate> > &candidates, bool recursive);
	void process(list_t &list, Candidate *candidate);
	void clear();
};

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleChain &chain);

#endif /* MPC_MODULE_CHAIN_H_ */

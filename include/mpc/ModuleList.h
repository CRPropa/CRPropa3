#ifndef MPC_MODULE_LIST_H_
#define MPC_MODULE_LIST_H_

#include <list>

#include "mpc/Candidate.h"
#include "mpc/Module.h"
#include "mpc/Source.h"

namespace mpc {

/**
 @class ModuleList
 @brief List of modules
 */
class ModuleList {
	std::list<ref_ptr<Module> > modules;

public:
	typedef std::list<ref_ptr<Module> >::iterator iterator;
	typedef std::list<ref_ptr<Module> >::const_iterator const_iterator;

	void add(Module *module);
	void process(Candidate *candidate);
	void run(Candidate *candidate, bool recursive);
	void run(Source *source, size_t count, bool recursive);

	const std::list<ref_ptr<Module> > &getModules() const;
};

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleList &list);

#endif /* MPC_MODULE_LIST_H_ */

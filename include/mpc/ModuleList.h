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
class ModuleList: public std::list<ref_ptr<Module> > {
public:
	typedef std::list<ref_ptr<Module> >::iterator iterator;
	typedef std::list<ref_ptr<Module> >::const_iterator const_iterator;

	void process(Candidate *candidate);
	void run(Candidate *candidate, bool recursive);
	void run(Source *source, size_t count, bool recursive);
};

} // namespace mpc

std::ostream &operator<<(std::ostream &out,
		const std::list<mpc::ref_ptr<mpc::Module> > &chain);

#endif /* MPC_MODULE_LIST_H_ */

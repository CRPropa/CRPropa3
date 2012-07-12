#ifndef MPC_MODULE_LIST_H_
#define MPC_MODULE_LIST_H_

#include "mpc/Candidate.h"
#include "mpc/Module.h"
#include "mpc/Source.h"
#include "mpc/AssocVector.h"
#include "mpc/Referenced.h"

#include <list>
#include <iostream>

namespace mpc {

/**
 @class ModuleList
 @brief List of modules
 */
class ModuleList: public Referenced {
public:
	typedef std::list<ref_ptr<Module> > module_list_t;
	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	ModuleList();
	virtual ~ModuleList();
	void setShowProgress(bool show);

	void add(Module* module);
	virtual void process(Candidate *candidate);
	void run(Candidate *candidate, bool recursive = true);
	void run(candidate_vector_t &candidates, bool recursive = true);
	void run(Source *source, size_t count, bool recursive = true);

	module_list_t &getModules();
	const module_list_t &getModules() const;

private:
	module_list_t modules;
	bool showProgress;
};

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleList &list);

#endif /* MPC_MODULE_LIST_H_ */

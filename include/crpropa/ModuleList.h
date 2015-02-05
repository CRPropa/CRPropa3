#ifndef CRPROPA_MODULE_LIST_H
#define CRPROPA_MODULE_LIST_H

#include "crpropa/Candidate.h"
#include "crpropa/Module.h"
#include "crpropa/Source.h"

#include <list>
#include <sstream>

namespace crpropa {

/**
 @class ModuleList
 @brief The simulation itself: A list of simulation modules
 */
class ModuleList: public Referenced {
public:
	typedef std::list<ref_ptr<Module> > module_list_t;
	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	ModuleList();
	virtual ~ModuleList();
	void setShowProgress(bool show = true); ///< activate a progress bar

	void add(Module* module);
	module_list_t &getModules();
	const module_list_t &getModules() const;

	/**
	 @class
	 @brief The simulation itself: A list of simulation modules
	 */
	void beginRun(); ///< call beginRun in all modules
	void endRun(); ///< call endRun in all modules
	void process(Candidate *candidate); ///< call process in all modules

	void run(Candidate *candidate, bool recursive = true);
	void run(candidate_vector_t &candidates, bool recursive = true);
	void run(SourceInterface *source, size_t count, bool recursive = true);

	std::string getDescription() const;
	void showModules() const;

private:
	module_list_t modules;
	bool showProgress;
};

} // namespace crpropa

#endif // CRPROPA_MODULE_LIST_H

#ifndef CRPROPA_MODULE_LIST_H
#define CRPROPA_MODULE_LIST_H

#include <algorithm>
#include <csignal>
#include <iostream>
#include <vector>
#include <exception>
#include <sstream>
#include <list>

#include "crpropa/Candidate.h"
#include "crpropa/Module.h"
#include "crpropa/Source.h"
#include "crpropa/module/Output.h"


namespace crpropa {

/**
 @class ModuleList
 @brief The simulation itself: A list of simulation modules
 */
class ModuleList: public Module {
public:
	typedef std::list<ref_ptr<Module> > module_list_t;
	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	ModuleList();
	virtual ~ModuleList();
	void setShowProgress(bool show = true); ///< activate a progress bar

	void add(ref_ptr<Module> module);
	void remove(std::size_t i);
	std::size_t size() const;
	ref_ptr<Module> operator[](const std::size_t i);

	void process(ref_ptr<Candidate> candidate) const; ///< call process in all modules

	void run(ref_ptr<Candidate> candidate, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a single candidate
	void run(ref_ptr<candidate_vector_t> candidates, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a candidate vector
	void run(ref_ptr<SourceInterface> source, size_t count, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a number of candidates from the given source	

	std::string getDescription() const;
	void showModules() const;
	
	/** iterator goodies */
	typedef module_list_t::iterator iterator;
	typedef module_list_t::const_iterator const_iterator;
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

	void setInterruptAction(ref_ptr<Output> action);
	void dumpCandidate(ref_ptr<Candidate> cand) const;

private:
	module_list_t modules;
	bool showProgress;
	ref_ptr<Output> interruptAction;
	bool haveInterruptAction = false;
	std::vector<int> notFinished; // list with not finished numbers of candidates
};

/**
 @class ModuleListRunner
 @brief Run the provided ModuleList when process is called.
 */
class ModuleListRunner: public Module {
private:
	ref_ptr<ModuleList> mlist;
public:

	ModuleListRunner(ref_ptr<ModuleList> mlist);
	void process(ref_ptr<Candidate> candidate) const; ///< call run of wrapped ModuleList
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_MODULE_LIST_H

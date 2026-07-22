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

 To use this module first create a instance over the constructor `ModuleList SIM();`.
 After that add some modules with `SIM.add(SomeModule)`
 When all module are added you can run the simulation over `SIM.run()`

 For a detailed description see the function descriptions.
 */
class ModuleList: public Module {
public:
	typedef std::list<ref_ptr<Module> > module_list_t;
	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	/// Constructor
	ModuleList();
	/// Destructor
	virtual ~ModuleList();
	/** Activates/Deactivates the progess bar
	 * If you want to see the progess of the simulation in stdout you can
	 * activate the progess bar with this function.
	 * @param show  Whether or not to show a progress bar 
	 */
	void setShowProgress(bool show = true);

	/** Function to add modules to the simulation
	 * With this function you can add any class derived from Module to the simulation.
	 * During ModuleList::run the Module::process function is called to process the simulated Candidate
	 * @param module  Module raw pointer
	 */
	void add(ref_ptr<Module> module);
	/** Removes module from index i in modules vector
	 * This function removes the module located at index i in the modules vector.
	 * The modules are added in order to the back of the modules vector.
	 * @param i  Index to remove 
	 */
	void remove(std::size_t i);
	/// @return Returns size of the modules vector, so the number of modules in ModuleList
	std::size_t size() const;
	/** Operator to get the module at index i
	 * @param i  Module index in modules vector
	 */
	ref_ptr<Module> operator[](const std::size_t i);

	/** process function
	 * The process function calls process in order on very module in the modules vector
	 * with the given Candidate reference pointer
	 * @param candidate  The candidate to call modules process on
	 */
	void process(ref_ptr<Candidate> candidate) const;

	/** Run function
	 * This run function does the full simulation for a single given candidate
	 * @param candidate  Candidate to use for the simulation
	 * @param recursive  Whether to also process possible created secondaries
	 * @param secondariesFirst  Whether to process the secondaries directly after the
	 * step where they are created or to only process them after the primary is fully finished
	 * @param waitForSecondaries  Whether to wait for the secondaries to be fully processed 
	 * before continuing with the primary candidate if secondariesFirst is set to true
	 */
	void run(ref_ptr<Candidate> candidate, bool recursive = true, bool secondariesFirst = true, bool waitForSecondaries = true);
	/** Run function
	 * This function calls the run function for a single Candidate on each entry of the given candidate_vector.
	 * If OpenMP parallelization is activated (default) the execution of the run function for each member
	 * of candidate_vector is parallelized. 
	 * @param candidates  Candidate vector
	 * @param recursive  Whether to also process possible created secondaries
	 * @param secondariesFirst  Whether to process the secondaries directly after the
	 * step where they are created or to only process them after the primary is fully finished
	 * @param waitForSecondaries  Whether to wait for the secondaries to be fully processed 
	 * before continuing with the primary candidate if secondariesFirst is set to true
	 */
	void run(ref_ptr<candidate_vector_t> candidates, bool recursive = true, bool secondariesFirst = false, bool waitForSecondaries = true);
	/** Run function
	 * This function generates count Candidates over the given source.
	 * Each generated Candidate is then handed over to ModuleList::run(ref_ptr<Candidate>)
	 * where it is simulated.
	 * If OpenMP parallelization is activated (default) the count primary Candidates are
	 * generated and simulated in parallel. 
	 * @param source  The source to generate Candidates from
	 * @param count  The number of Candidates to simulate
	 * @param recursive  Whether to also process possible created secondaries
	 * @param secondariesFirst  Whether to process the secondaries directly after the
	 * step where they are created or to only process them after the primary is fully finished
	 * @param waitForSecondaries  Whether to wait for the secondaries to be fully processed 
	 * before continuing with the primary candidate if secondariesFirst is set to true
	 */
	void run(ref_ptr<SourceInterface> source, size_t count, bool recursive = true, bool secondariesFirst = false, bool waitForSecondaries = true);

	/// @return Returns the description string for all modules
	std::string getDescription() const;
	/// @brief Prints the getDescription string to stdout
	void showModules() const;
	
	/* iterator goodies */
	typedef module_list_t::iterator iterator;
	typedef module_list_t::const_iterator const_iterator;
	iterator begin(); /**< Begin iterator of modules vector */
	const_iterator begin() const; /**< Constant begin iterator of modules vector */
	iterator end(); /**< End iterator of modules vector */
	const_iterator end() const; /**< Constant end iterator of modules vector */

	/** Sets outputmodule that should be called on interruption
	 * If CRPropa receives an interruption signal (4) the here set Output module is called
	 * @param action  The Output module to call when CRPropa is interrupted
	 */
	void setInterruptAction(ref_ptr<Output> action);
	/** Function to dump Candidates manually
	 * With this function it is possible to dump the currently still active simulated Candidates
	 * manually. When called the with ModuleList::setInterruptAction set module is called on
	 * all still active Candidates.
	 * @param Candidate  Candidate to dump 
	 */
	void dumpCandidate(ref_ptr<Candidate> cand) const;

private:
	module_list_t modules; /**< Vector containing reference pointer of all modules added with ModuleList::add */
	bool showProgress; /**< Whether to show the progess bar or not */
	ref_ptr<Output> interruptAction; /**< The output module to call when CRPropa receives a userinterrupt signal */
	bool haveInterruptAction = false; /**< Whether a interrupt action is set or not */
	std::vector<int> notFinished; /**< List with not finished numbers of candidates */
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

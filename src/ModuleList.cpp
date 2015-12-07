#include "crpropa/ModuleList.h"
#include "crpropa/ProgressBar.h"

#if _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <signal.h>
#ifndef sighandler_t
typedef void (*sighandler_t)(int);
#endif

using namespace std;

namespace crpropa {

bool g_cancel_signal_flag = false;
void g_cancel_signal_callback(int sig) {
	std::cerr << "crpropa::ModuleList: SIGINT/SIGTERM received" << std::endl;
	g_cancel_signal_flag = true;
}

ModuleList::ModuleList() :
		showProgress(false) {
}

ModuleList::~ModuleList() {
}

void ModuleList::setShowProgress(bool show) {
	showProgress = show;
}

void ModuleList::add(Module *module) {
	modules.push_back(module);
}

void ModuleList::beginRun() {
	module_list_t::iterator m;
	for (m = modules.begin(); m != modules.end(); m++)
		(*m)->beginRun();
}

void ModuleList::process(Candidate *candidate) const {
	module_list_t::const_iterator m;
	for (m = modules.begin(); m != modules.end(); m++)
		(*m)->process(candidate);
}

void ModuleList::endRun() {
	module_list_t::iterator m;
	for (m = modules.begin(); m != modules.end(); m++)
		(*m)->endRun();
}

void ModuleList::run(Candidate *candidate, bool recursive) {
	while (candidate->isActive() && !g_cancel_signal_flag)
		process(candidate);

	// propagate secondaries
	if (recursive) {
		for (size_t i = 0; i < candidate->secondaries.size(); i++) {
			if (g_cancel_signal_flag)
				break;
			run(candidate->secondaries[i], recursive);
		}
	}


}

void ModuleList::run(candidate_vector_t &candidates, bool recursive) {
	size_t count = candidates.size();

#if _OPENMP
	std::cout << "crpropa::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar progressbar(count);

	if (showProgress) {
		progressbar.start("Run ModuleList");
	}

	g_cancel_signal_flag = false;
	sighandler_t old_sigint_handler = ::signal(SIGINT,
			g_cancel_signal_callback);
	sighandler_t old_sigterm_handler = ::signal(SIGTERM,
			g_cancel_signal_callback);

	beginRun();  // call beginRun in all modules

#pragma omp parallel for schedule(static, 1000)
	for (size_t i = 0; i < count; i++) {
		if (g_cancel_signal_flag)
			continue;

		try {
			run(candidates[i], recursive);
		} catch (std::exception &e) {
			std::cerr << "Exception in crpropa::ModuleList::run: " << std::endl;
			std::cerr << e.what() << std::endl;
		}

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	endRun();  // call endRun in all modules

	::signal(SIGINT, old_sigint_handler);
	::signal(SIGTERM, old_sigterm_handler);
}

void ModuleList::run(SourceInterface *source, size_t count, bool recursive) {

#if _OPENMP
	std::cout << "crpropa::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar progressbar(count);

	if (showProgress) {
		progressbar.start("Run ModuleList");
	}

	g_cancel_signal_flag = false;
	sighandler_t old_signal_handler = ::signal(SIGINT,
			g_cancel_signal_callback);

	beginRun();  // call beginRun in all modules

#pragma omp parallel for schedule(static, 1000)
	for (size_t i = 0; i < count; i++) {
		if (g_cancel_signal_flag)
			continue;

		ref_ptr<Candidate> candidate;
		
		try {
			candidate = source->getCandidate();
		} catch (std::exception &e) {
			std::cerr << "Exception in crpropa::ModuleList::run: source->getCandidate" << std::endl;
			std::cerr << e.what() << std::endl;
		}

		if (candidate.valid()) {
			try {
				run(candidate, recursive);
			} catch (std::exception &e) {
				std::cerr << "Exception in crpropa::ModuleList::run: " << std::endl;
				std::cerr << e.what() << std::endl;
			}
		}

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	endRun();  // call endRun in all modules

	::signal(SIGINT, old_signal_handler);
}

ModuleList::module_list_t &ModuleList::getModules() {
	return modules;
}

const ModuleList::module_list_t &ModuleList::getModules() const {
	return modules;
}

std::string ModuleList::getDescription() const {
	std::stringstream ss;
	ss << "ModuleList\n";
	crpropa::ModuleList::module_list_t::const_iterator m;
	for (m = modules.begin(); m != modules.end(); m++)
		ss << "  " << (*m)->getDescription() << "\n";
	return ss.str();
}

void ModuleList::showModules() const {
	std::cout << getDescription();
}

ModuleListRunner::ModuleListRunner(ModuleList *mlist) : mlist(mlist) {
}

void ModuleListRunner::process(Candidate *candidate) const {
	if (mlist.valid())
		mlist->run(candidate);
}

std::string ModuleListRunner::getDescription() const {
	std::stringstream ss;
	ss << "ModuleListRunner\n";
	if (mlist.valid())
		ss << mlist->getDescription();
	return ss.str();
};

} // namespace crpropa

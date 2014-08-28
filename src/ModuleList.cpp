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

void ModuleList::process(Candidate *candidate) {
	module_list_t::iterator iEntry = modules.begin();
	while (iEntry != modules.end()) {
		ref_ptr<Module> &module = *iEntry;
		iEntry++;
		module->process(candidate);
	}
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
	sighandler_t old_signal_handler = ::signal(SIGINT,
			g_cancel_signal_callback);

#pragma omp parallel for schedule(static, 1000)
	for (size_t i = 0; i < count; i++) {
		if (g_cancel_signal_flag)
			continue;

		run(candidates[i], recursive);

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	::signal(SIGINT, old_signal_handler);
}

void ModuleList::run(Source *source, size_t count, bool recursive) {

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

#pragma omp parallel for schedule(static, 1000)
	for (size_t i = 0; i < count; i++) {
		if (g_cancel_signal_flag)
			continue;

		ref_ptr<Candidate> candidate = source->getCandidate();
		run(candidate, recursive);

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

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
	crpropa::ModuleList::module_list_t::const_iterator it;
	for (it = modules.begin(); it != modules.end(); ++it) {
		const crpropa::ref_ptr<crpropa::Module> &m = *it;
		ss << "  " << m->getDescription() << "\n";
	}
	return ss.str();
}

void ModuleList::showModules() const {
	std::cout << getDescription();
}

} // namespace crpropa

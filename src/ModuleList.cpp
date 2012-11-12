#include "mpc/ModuleList.h"
#include "mpc/ProgressBar.h"

#include <omp.h>
#include <signal.h>
#include <algorithm>

using namespace std;

namespace mpc {

bool g_cancel_signal_flag = false;
sighandler_t g_cancel_signal_backup = NULL;
void g_cancel_signal_callback(int sig) {
	g_cancel_signal_flag = true;
	::signal(SIGINT, SIG_DFL );
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
	if (recursive)
		for (size_t i = 0; i < candidate->secondaries.size(); i++)
			run(candidate->secondaries[i], recursive);
}

void ModuleList::run(candidate_vector_t &candidates, bool recursive) {
	size_t count = candidates.size();

#if _OPENMP
	std::cout << "mpc::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar progressbar(count);

	if (showProgress) {
		progressbar.start("Run ModuleList");
	}

	g_cancel_signal_flag = false;
	g_cancel_signal_backup = ::signal(SIGINT, g_cancel_signal_callback);

#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
		if (!g_cancel_signal_flag)
			run(candidates[i], recursive);

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	::signal(SIGINT, g_cancel_signal_backup);
}

void ModuleList::run(Source *source, size_t count, bool recursive) {

#if _OPENMP
	std::cout << "mpc::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar progressbar(count);

	if (showProgress) {
		progressbar.start("Run ModuleList");
	}

	g_cancel_signal_flag = false;
	::signal(SIGINT, g_cancel_signal_callback);

#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
		ref_ptr<Candidate> candidate = source->getCandidate();
		if (!g_cancel_signal_flag)
			run(candidate, recursive);

		if (showProgress)
#pragma omp critical(progressbarUpdate)
			progressbar.update();
	}

	::signal(SIGINT, g_cancel_signal_backup);
}

ModuleList::module_list_t &ModuleList::getModules() {
	return modules;
}

const ModuleList::module_list_t &ModuleList::getModules() const {
	return modules;
}

void ModuleList::showModules() const {
	mpc::ModuleList::module_list_t::const_iterator iEntry;

	iEntry = getModules().begin();
	while (iEntry != getModules().end()) {
		const mpc::ref_ptr<mpc::Module> &entry = *iEntry;
		iEntry++;
		std::cout << "  " << entry->getDescription() << "\n";
	}
}

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleList &list) {
	mpc::ModuleList::module_list_t::const_iterator iEntry;

	iEntry = list.getModules().begin();
	while (iEntry != list.getModules().end()) {
		const mpc::ref_ptr<mpc::Module> &entry = *iEntry;
		iEntry++;
		out << "  " << entry->getDescription() << "\n";
	}
	return out;
}

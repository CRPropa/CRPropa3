#include "mpc/ModuleList.h"
#include "mpc/ProgressBar.h"

#include <omp.h>
#include <algorithm>

using namespace std;

namespace mpc {

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
	while (candidate->isActive()) {
		process(candidate);
	}

	// propagate secondaries
	if (recursive) {
		for (size_t i = 0; i < candidate->secondaries.size(); i++)
			run(candidate->secondaries[i], recursive);
	}

}

void ModuleList::run(candidate_vector_t &candidates, bool recursive) {
	size_t count = candidates.size();
	size_t pc = std::max(count / 100, size_t(1));
#if _OPENMP
	std::cout << "mpc::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar *progressbar = NULL;
	if (showProgress) {
		progressbar = new ProgressBar("Run ModuleList", count);
	}

#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
		//if ((showProgress) && (i % pc == 0))
		//	std::cout << i / pc << "% - " << i << std::endl;
		run(candidates[i], recursive);

		if (progressbar)
#pragma omp critical(progressbarUpdate)
			progressbar->update();
	}
}

void ModuleList::run(Source *source, size_t count, bool recursive) {
	size_t pc = std::max(count / 100, size_t(1));

#if _OPENMP
	std::cout << "mpc::ModuleList: Number of Threads: " << omp_get_max_threads() << std::endl;
#endif

	ProgressBar *progressbar = NULL;
	if (showProgress) {
		progressbar = new ProgressBar("Run ModuleList", count);
	}

#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
		ParticleState state;
		source->prepare(state);
		ref_ptr<Candidate> candidate = new Candidate(state);
		run(candidate, recursive);
		if (progressbar)
#pragma omp critical(progressbarUpdate)
			progressbar->update();
	}
}

ModuleList::module_list_t &ModuleList::getModules() {
	return modules;
}

const ModuleList::module_list_t &ModuleList::getModules() const {
	return modules;
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

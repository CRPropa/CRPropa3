#include "mpc/ModuleList.h"

#include <omp.h>

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
	size_t cent = count / 100;
	if (cent == 0)
		cent = 1;
	size_t pc = 0;
#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
#if _OPENMP
		if (i == 0) {
			std::cout << "Number of Threads: " << omp_get_num_threads()
			<< std::endl;
		}
#endif
		if (showProgress && (i % cent == 0)) {
			std::cout << pc << "% - " << i << std::endl;
			pc++;
		}
		run(candidates[i], recursive);
	}
}

void ModuleList::run(Source *source, size_t count, bool recursive) {
	size_t cent = count / 100;
	if (cent == 0)
		cent = 1;
	size_t pc = 0;
#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
#if _OPENMP
		if (i == 0) {
			std::cout << "Number of Threads: " << omp_get_num_threads()
			<< std::endl;
		}
#endif
		if (showProgress && (i % cent == 0)) {
			std::cout << pc << "% - " << i << std::endl;
			pc++;
		}
		ParticleState state;
		source->prepare(state);
		ref_ptr<Candidate> candidate = new Candidate(state);
		run(candidate, recursive);
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

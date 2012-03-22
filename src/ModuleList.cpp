#include "mpc/ModuleList.h"

using namespace std;

namespace mpc {

void ModuleList::add(Module *module) {
	modules.push_back(module);
}

void ModuleList::process(Candidate *candidate) {
	iterator iEntry = modules.begin();
	while (iEntry != modules.end()) {
		ref_ptr<Module> &module = *iEntry;
		iEntry++;
		module->process(candidate);
	}
}

void ModuleList::run(Candidate *candidate, bool recursive) {
	while (candidate->getStatus() == Candidate::Active) {
		process(candidate);
	}

	// propagate secondaries
	if (recursive) {
		for (size_t i = 0; i < candidate->secondaries.size(); i++)
			run(candidate->secondaries[i], recursive);
	}

}

void ModuleList::run(Source *source, size_t count, bool recursive) {
#pragma omp parallel for schedule(dynamic, 1000)
	for (size_t i = 0; i < count; i++) {
		ParticleState state;
		source->prepare(state);
		ref_ptr<Candidate> candidate = new Candidate(state);
		run(candidate, recursive);
	}
}

const std::list<ref_ptr<Module> > &ModuleList::getModules() const {
	return modules;
}

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleList &list) {
	mpc::ModuleList::const_iterator iEntry;

	iEntry = list.getModules().begin();
	while (iEntry != list.getModules().end()) {
		const mpc::ref_ptr<mpc::Module> &entry = *iEntry;
		iEntry++;
		out << "  " << entry->getDescription() << "\n";
	}
	return out;
}

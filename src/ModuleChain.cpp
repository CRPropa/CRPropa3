#include "mpc/ModuleChain.h"

namespace mpc {

const ModuleChain::list_t &ModuleChain::getStartModules() const {
	return startModules;
}
const ModuleChain::list_t &ModuleChain::getMainModules() const {
	return mainModules;
}
const ModuleChain::list_t &ModuleChain::getEndModules() const {
	return endModules;
}

void ModuleChain::add(size_t priority, Module *module) {
	if (module == 0)
		return;

	list_entry_t entry;
	entry.first = priority;
	entry.second = module;

	if (priority == Priority::Start) {
		startModules.push_back(entry);
	} else if (priority == Priority::End) {
		endModules.push_back(entry);
	} else {
		mainModules.push_back(entry);
		mainModules.sort();
	}
}

void ModuleChain::clear() {
	startModules.clear();
	mainModules.clear();
	endModules.clear();
}

void ModuleChain::process(list_t &list, Candidate *candidate) {
	list_t::iterator iEntry = list.begin();
	while (iEntry != list.end()) {
		list_entry_t &entry = *iEntry;
		iEntry++;
		entry.second->process(candidate);
	}
}

void ModuleChain::process(Candidate *candidate) {
	if (mainModules.size() == 0)
		return;

	process(startModules, candidate);

	while (candidate->getStatus() == Candidate::Active) {
		if (mainModules.size() == 0)
			break;
		process(mainModules, candidate);
	}

	process(endModules, candidate);
}

void ModuleChain::process(std::vector<ref_ptr<Candidate> > &candidates,
		bool recursive) {
#pragma omp parallel for
	for (size_t i = 0; i < candidates.size(); i++) {
		Candidate *candidate = candidates[i];
		if (candidate->getStatus() != Candidate::Active)
			continue;
		process(candidate);

		// propagate secondaries
		if (recursive)
			process(candidate->secondaries, recursive);
	}
}

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleChain &chain) {
	mpc::ModuleChain::list_t::const_iterator iEntry;

	out << "Start Modules" << "\n";
	iEntry = chain.getStartModules().begin();
	while (iEntry != chain.getStartModules().end()) {
		const mpc::ModuleChain::list_entry_t &entry = *iEntry;
		iEntry++;
		out << "  " << entry.first << " -> " << entry.second->getDescription()
				<< "\n";
	}

	out << "\nMain Modules" << "\n";
	iEntry = chain.getMainModules().begin();
	while (iEntry != chain.getMainModules().end()) {
		const mpc::ModuleChain::list_entry_t &entry = *iEntry;
		iEntry++;
		out << "  " << entry.first << " -> " << entry.second->getDescription()
				<< "\n";
	}

	out << "\nEnd Modules" << "\n";
	iEntry = chain.getEndModules().begin();
	while (iEntry != chain.getEndModules().end()) {
		const mpc::ModuleChain::list_entry_t &entry = *iEntry;
		iEntry++;
		out << "  " << entry.first << " -> " << entry.second->getDescription()
				<< "\n";
	}
	return out;
}

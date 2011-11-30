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

	check();
}

void ModuleChain::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {
	list_t::iterator iStartEntry = startModules.begin();
	while (iStartEntry != startModules.end()) {
		list_entry_t &entry = *iStartEntry;
		iStartEntry++;

		entry.second->process(candidate, secondaries);
	}

	while (candidate->getStatus() == Candidate::Active) {
		list_t::iterator iEntry = mainModules.begin();
		while (iEntry != mainModules.end()) {
			list_entry_t &entry = *iEntry;
			iEntry++;

			entry.second->process(candidate, secondaries);
		}
	}

	list_t::iterator iEndEntry = endModules.begin();
	while (iEndEntry != endModules.end()) {
		list_entry_t &entry = *iEndEntry;
		iEndEntry++;

		entry.second->process(candidate, secondaries);
	}
}

void ModuleChain::clear() {
	startModules.clear();
	mainModules.clear();
	endModules.clear();
}

void ModuleChain::check() {
	size_t integratorCount = 0;

	list_t::iterator iEntry = mainModules.begin();
	while (iEntry != mainModules.end()) {
		list_entry_t &entry = *iEntry;
		iEntry++;

		if (entry.first == Priority::Integration) {
			integratorCount++;
		}
	}

	if (integratorCount > 1) {
		std::cerr << "Warning: more than one integration feature present."
				<< std::endl;
	}
}

void ModuleChain::process(std::vector<Candidate *> &candidates) {
	bool haveActive = true;
	while (haveActive) {
		haveActive = false;
		for (size_t i = 0; i < candidates.size(); i++) {
			Candidate *candidate = candidates[i];
			if (candidate->getStatus() == Candidate::Active) {
				process(candidate, candidates);
				haveActive = true;
			}
		}
	}
}

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::ModuleChain &chain) {
	mpc::ModuleChain::list_t::const_iterator iEntry =
			chain.getMainModules().begin();
	while (iEntry != chain.getMainModules().end()) {
		const mpc::ModuleChain::list_entry_t &entry = *iEntry;
		iEntry++;

		out << entry.first << " -> " << entry.second->getDescription() << "\n";
	}
	return out;
}

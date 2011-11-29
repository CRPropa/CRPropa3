#ifndef MODULE_CHAIN_H_
#define MODULE_CHAIN_H_

#include <list>
#include <typeinfo>
#include <assert.h>

#include "mpc/SharedPointer.h"
#include "mpc/Candidate.h"
#include "mpc/Module.h"

namespace mpc {

struct Priority {
	enum Enum {
		Start = 0,
		BeforeIntegration = 20,
		Integration = 25,
		AfterIntegration = 30,
		BeforeInteraction = 45,
		Interaction = 50,
		AfterInteraction = 55,
		BeforeCommit = 70,
		Commit = 75,
		AfterCommit = 80,
		End = 100
	};
};

class ModuleChain {
	struct ModuleEntry {
		size_t priority;
		shared_ptr<Module> module;
		bool operator<(ModuleEntry const& rhs) const {
			return (priority < rhs.priority);
		}
	};
	std::list<ModuleEntry> startModules;
	std::list<ModuleEntry> mainModules;
	std::list<ModuleEntry> endModules;

	void check() {
		size_t integratorCount = 0;

		std::list<ModuleEntry>::iterator iEntry = mainModules.begin();
		while (iEntry != mainModules.end()) {
			ModuleEntry &entry = *iEntry;
			iEntry++;

			if (entry.priority == Priority::Integration) {
				integratorCount++;
			}
		}

		if (integratorCount > 1) {
			std::cerr << "Warning: more than one integration feature present."
					<< std::endl;
		}
	}
public:

	void add(size_t priority, shared_ptr<Module> feature) {
		ModuleEntry entry;
		entry.priority = priority;
		entry.module = feature;

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

	void apply(Candidate &candidate) {
		std::list<ModuleEntry>::iterator iStartEntry = startModules.begin();
		while (iStartEntry != startModules.end()) {
			ModuleEntry &entry = *iStartEntry;
			iStartEntry++;

			entry.module->apply(candidate);
		}

		while (candidate.getStatus() == Candidate::Active) {
			std::list<ModuleEntry>::iterator iEntry = mainModules.begin();
			while (iEntry != mainModules.end()) {
				ModuleEntry &entry = *iEntry;
				iEntry++;

				entry.module->apply(candidate);
			}
		}

		std::list<ModuleEntry>::iterator iEndEntry = endModules.begin();
		while (iEndEntry != endModules.end()) {
			ModuleEntry &entry = *iEndEntry;
			iEndEntry++;

			entry.module->apply(candidate);
		}
	}

	void print(std::ostream &out = std::cout) {
		std::list<ModuleEntry>::iterator iEntry = mainModules.begin();
		while (iEntry != mainModules.end()) {
			ModuleEntry &entry = *iEntry;
			iEntry++;

			out << entry.priority << " -> " << entry.module->getDescription()
					<< "\n";
		}
		out.flush();
	}

	void clear() {
		startModules.clear();
		mainModules.clear();
		endModules.clear();
	}
};

} // namespace mpc

#endif /* MODULE_CHAIN_H_ */

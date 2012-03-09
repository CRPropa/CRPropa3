#include "mpc/ModuleList.h"

using namespace std;

namespace mpc {

void ModuleList::process(Candidate *candidate) {
	iterator iEntry = begin();
	while (iEntry != end()) {
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

} // namespace mpc

std::ostream &operator<<(std::ostream &out,
		const std::list<mpc::ref_ptr<mpc::Module> > &modules) {
	std::list<mpc::ref_ptr<mpc::Module> >::const_iterator iEntry;

	iEntry = modules.begin();
	while (iEntry != modules.end()) {
		const mpc::ref_ptr<mpc::Module> &entry = *iEntry;
		iEntry++;
		out << "  " << entry->getDescription() << "\n";
	}
	return out;
}

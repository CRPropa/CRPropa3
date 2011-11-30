#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "mpc/Module.h"

namespace mpc {

class CandidateOutput: public Module {
	std::string getDescription() const {
		return "CandidateOutput";
	}

	void apply(Candidate &candidate) {
		std::cout << "Age: " << candidate.getTrajectoryLength() / Mpc << std::endl;
		std::cout << "  CurrentStep: " << candidate.getCurrentStep() * c_light / Mpc << std::endl;
		std::cout << "  NextStep:    " << candidate.getNextStep() * c_light / Mpc << std::endl;
	}
};

} // namespace mpc

#endif /* OUTPUT_H_ */

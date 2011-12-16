#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "mpc/Module.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace mpc {

class CandidateOutput: public Module {
public:
	std::string getDescription() const {
		return "CandidateOutput";
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		std::cout << "Age: " << candidate->getTrajectoryLength() / Mpc << std::endl;
		std::cout << "  CurrentStep: " << candidate->getCurrentStep() * c_light / Mpc << std::endl;
		std::cout << "  NextStep:    " << candidate->getNextStep() * c_light / Mpc << std::endl;
	}
};

class TrajectoryOutput: public Module {
private:
	std::ofstream outfile;

public:
	std::string getDescription() const {
		return "TrajectoryOutput";
	}

	TrajectoryOutput(std::string outFileName) {
		outfile.open(outFileName.c_str());
		outfile << "Age [Mpc], HepId, E [EeV], posX, posY, posZ, dirX, dirY, dirZ\n";
	}

	~TrajectoryOutput() {
		outfile.close();
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		outfile << candidate->getTrajectoryLength() / Mpc << ", ";
		outfile << candidate->current.getId() << ", ";
		outfile << candidate->current.getEnergy() / EeV << ", ";
		Vector3 pos = candidate->current.getPosition() / Mpc;
		outfile << pos.x() << ", ";
		outfile << pos.y() << ", ";
		outfile << pos.z() << ", ";
		Vector3 dir = candidate->current.getDirection();
		outfile << dir.x() << ", ";
		outfile << dir.y() << ", ";
		outfile << dir.z() << "\n";
	}
};

} // namespace mpc

#endif /* OUTPUT_H_ */

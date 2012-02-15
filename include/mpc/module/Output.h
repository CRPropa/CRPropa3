#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "mpc/Module.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace mpc {

std::string getOutputString(ParticleState particle) {
	std::stringstream ss;
	ss << particle.getId() << ", ";
	ss << particle.getEnergy() / EeV << ", ";
	Vector3 pos = particle.getPosition() / Mpc;
	ss << pos.x() << ", ";
	ss << pos.y() << ", ";
	ss << pos.z() << ", ";
	Vector3 dir = particle.getDirection();
	ss << dir.x() << ", ";
	ss << dir.y() << ", ";
	ss << dir.z();
	return ss.str();
}

/**
 @class TrajectoryOutput
 @brief Saves trajectories to CSV file.
 */
class TrajectoryOutput: public Module {
private:
	std::ofstream outfile;

public:
	TrajectoryOutput(std::string name) {
		outfile.open(name.c_str());
		outfile << "# Age, HepId, E, posX, posY, posZ, dirX, dirY, dirZ\n";
	}

	~TrajectoryOutput() {
		outfile.close();
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		outfile << candidate->getTrajectoryLength() / Mpc << ", "
				<< getOutputString(candidate->current);
	}

	std::string getDescription() const {
		return "TrajectoryOutput";
	}
};

/**
 @class FinalOutput
 @brief Saves particles to a CSV file.
 */
class FinalOutput: public Module {
private:
	std::ofstream outfile;
	Candidate::Status flag;

public:
	FinalOutput(std::string name, Candidate::Status flag) {
		this->flag = flag;
		outfile.open(name.c_str());
		outfile << "(initial) Age, Id, E, x, y, z, dirX, dirY, dirZ, ";
		outfile << "(final) Age, Id, E, x, y, z, dirX, dirY, dirZ\n";
	}

	~FinalOutput() {
		outfile.close();
	}

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		std::cout << "[FinalOutput] process" << std::endl;
		if (candidate->getStatus() != flag)
			return;
		// initial state
		outfile << "0, ";
		outfile << getOutputString(candidate->initial);
		outfile << ", ";
		// current state
		outfile << candidate->getTrajectoryLength() / Mpc << ", ";
		outfile << getOutputString(candidate->current);
		outfile << "\n";
		outfile.close();
	}

	std::string getDescription() const {
		return "FinishedOutput";
	}
};

/**
 @class ShellOutput
 @brief Output of the candidate to the shell.
 */
class ShellOutput: public Module {
public:
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries) {
		std::cout << std::fixed << std::showpoint << std::setprecision(2)
				<< std::setw(6);
		std::cout << candidate->getTrajectoryLength() / Mpc << " Mpc,  ";
		std::cout << candidate->current.getId() << ",  ";
		std::cout << candidate->current.getEnergy() / EeV << " EeV,  ";
		std::cout << candidate->current.getPosition() / Mpc << " Mpc";
		std::cout << std::endl;
	}

	std::string getDescription() const {
		return "ShellOutput";
	}
};

} // namespace mpc

#endif /* OUTPUT_H_ */

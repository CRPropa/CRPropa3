#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "mpc/Module.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace mpc {

std::string getOutputString(ParticleState particle) {
	std::stringstream s;
	s << particle.getId() << ", ";
	s << particle.getEnergy() / EeV << ", ";
	Vector3 pos = particle.getPosition() / Mpc;
	s << pos.x() << ", ";
	s << pos.y() << ", ";
	s << pos.z() << ", ";
	Vector3 dir = particle.getDirection();
	s << dir.x() << ", ";
	s << dir.y() << ", ";
	s << dir.z();
	return s.str();
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
		outfile << "# Age, HepId, E, posX, posY, posZ, dirX, dirY, dirZ, event\n";
	}

	~TrajectoryOutput() {
		outfile.close();
	}

	void process(Candidate *candidate) {
		outfile << candidate->getTrajectoryLength() / Mpc << ", "
				<< getOutputString(candidate->current) << "\n";
	}

	std::string getDescription() const {
		return "Trajectory output";
	}
};

/**
 @class FlagOutput
 @brief Saves particles with a given flag to a CSV file.
 */
class FlagOutput: public Module {
private:
	std::ofstream outfile;
	Candidate::Status flag;

public:
	FlagOutput(std::string name, Candidate::Status flag) {
		this->flag = flag;
		outfile.open(name.c_str());
		outfile << "(initial) Age, Id, E, x, y, z, dirX, dirY, dirZ, ";
		outfile << "(flagged) Age, Id, E, x, y, z, dirX, dirY, dirZ\n";
	}

	~FlagOutput() {
		outfile.close();
	}

	void process(Candidate *candidate) {
		if (candidate->getStatus() != flag)
			return;
		// initial state
		outfile << "0, ";
		outfile << getOutputString(candidate->initial);
		outfile << ", ";
		// final state
		outfile << candidate->getTrajectoryLength() / Mpc << ", ";
		outfile << getOutputString(candidate->current);
		outfile << "\n";
	}

	std::string getDescription() const {
		return "Output, flag: ";
	}
};

/**
 @class ShellOutput
 @brief Write properties of the candidate to the shell.
 */
class ShellOutput: public Module {
public:
	void process(Candidate *candidate) {
		std::cout << std::fixed << std::showpoint << std::setprecision(2)
				<< std::setw(6);
		std::cout << candidate->getTrajectoryLength() / Mpc << " Mpc,  ";
		std::cout << candidate->current.getId() << ",  ";
		std::cout << candidate->current.getEnergy() / EeV << " EeV,  ";
		std::cout << candidate->current.getPosition() / Mpc << " Mpc, Status: ";
		std::cout << candidate->getStatus();
		std::cout << std::endl;
	}

	std::string getDescription() const {
		return "ShellOutput";
	}
};

} // namespace mpc

#endif /* OUTPUT_H_ */

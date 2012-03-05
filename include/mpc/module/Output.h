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
 @class FinalOutput
 @brief Saves particles with given flag to a CSV file.
 */
class FlaggedOutput: public Module {
private:
	std::ofstream outfile;
	Candidate::Status flag;

public:
	FlaggedOutput(std::string name, Candidate::Status flag) {
		this->flag = flag;
		outfile.open(name.c_str());
		outfile << "(initial) Age, Id, E, x, y, z, dirX, dirY, dirZ, ";
		outfile << "(flagged) Age, Id, E, x, y, z, dirX, dirY, dirZ\n";
	}

	~FlaggedOutput() {
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

class ChargeMassEngergyOutput: public Module {
private:
	Candidate::Status flag;
	std::ofstream outfile;

public:
	ChargeMassEngergyOutput(std::string filename, Candidate::Status flag) {
		this->flag = flag;
		outfile.open(filename.c_str());
		outfile << "Z, A, E" << std::endl;
	}

	~ChargeMassEngergyOutput() {
		outfile.close();
	}

	void process(Candidate *candidate) {
		if (candidate->getStatus() != flag)
			return;
		double Z = candidate->current.getChargeNumber();
		double A = candidate->current.getMassNumber();
		double E = candidate->current.getEnergy();
		outfile << Z << ", " << A << ", " << E << std::endl;
	}

	std::string getDescription() const {
		return "ChargeMassEnergyOutput";
	}
};

/**
 @class ShellOutput
 @brief Output of the candidate to the shell.
 */
class ShellOutput: public Module {
public:
	void process(Candidate *candidate) {
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

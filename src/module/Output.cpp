#include "mpc/module/Output.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <kiss/convert.h>

namespace mpc {

TrajectoryOutput::TrajectoryOutput(std::string name) {
	outfile.open(name.c_str());
	outfile << "# Age, HepId, E, posX, posY, posZ, dirX, dirY, dirZ, event\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	outfile.close();
}

void TrajectoryOutput::process(Candidate *candidate) const {
#pragma omp critical
	{
		outfile << candidate->getTrajectoryLength() / Mpc << ", ";
		outfile << candidate->current.getId() << ", ";
		outfile << candidate->current.getEnergy() / EeV << ", ";
		outfile << candidate->current.getPosition().x() / Mpc << ", ";
		outfile << candidate->current.getPosition().y() / Mpc << ", ";
		outfile << candidate->current.getPosition().z() / Mpc << ", ";
		outfile << candidate->current.getDirection().x() << ", ";
		outfile << candidate->current.getDirection().y() << ", ";
		outfile << candidate->current.getDirection().y() << ", ";
		outfile << std::endl;
	}
}

std::string TrajectoryOutput::getDescription() const {
	return "Trajectory output";
}

ConditionalOutput::ConditionalOutput(std::string filename,
		std::string propName) {
	propertyName = propName;
	outfile.open(filename.c_str());
	outfile
			<< "id, x, y, z, E, phi, theta, distance, i_id, i_x, i_y, i_z, i_E, i_phi, i_theta"
			<< std::endl;
}

ConditionalOutput::~ConditionalOutput() {
	outfile.close();
}

void ConditionalOutput::process(Candidate *candidate) const {
	if (candidate->hasProperty(propertyName)) {
#pragma omp critical
		{
			outfile << candidate->current.getId() << ", ";
			outfile << candidate->current.getPosition().x() / Mpc << ", ";
			outfile << candidate->current.getPosition().y() / Mpc << ", ";
			outfile << candidate->current.getPosition().z() / Mpc << ", ";
			outfile << candidate->current.getEnergy() / EeV << ", ";
			outfile << candidate->current.getDirection().phi() << ", ";
			outfile << candidate->current.getDirection().theta() << ", ";
			outfile << candidate->getTrajectoryLength() / Mpc << ", ";
			outfile << candidate->initial.getId() << ", ";
			outfile << candidate->initial.getPosition().x() / Mpc << ", ";
			outfile << candidate->initial.getPosition().y() / Mpc << ", ";
			outfile << candidate->initial.getPosition().z() / Mpc << ", ";
			outfile << candidate->initial.getEnergy() / EeV << ", ";
			outfile << candidate->initial.getDirection().phi() << ", ";
			outfile << candidate->initial.getDirection().theta();
			outfile << std::endl;
		}
	}
	candidate->removeProperty(propertyName);
}

std::string ConditionalOutput::getDescription() const {
	return "ConditionalOutput, condition: " + propertyName;
}

void ShellOutput::process(Candidate *candidate) const {
#pragma omp critical
	{
		std::cout << std::fixed << std::showpoint << std::setprecision(2)
				<< std::setw(6);
		std::cout << candidate->getTrajectoryLength() / Mpc << " Mpc,  ";
		std::cout << candidate->current.getId() << ",  ";
		std::cout << candidate->current.getEnergy() / EeV << " EeV,  ";
		std::cout << candidate->current.getPosition() / Mpc << " Mpc, Status: ";
		std::cout << std::endl;
	}
}

std::string ShellOutput::getDescription() const {
	return "ShellOutput";
}

} // namespace mpc

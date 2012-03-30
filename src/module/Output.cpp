#include "mpc/module/Output.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>

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
	char buffer[128];
	size_t pos = 0;
	pos += ::sprintf(buffer + pos, "%f, %d",
			candidate->getTrajectoryLength() / Mpc, candidate->current.getId());
	Vector3 position = candidate->current.getPosition() / Mpc;
	pos += ::sprintf(buffer + pos, ", %f, %f, %f", position.x(), position.y(),
			position.z());
	const Vector3 &dir = candidate->current.getDirection();
	pos += ::sprintf(buffer + pos, ", %f, %f, %f\n", dir.x(), dir.y(), dir.z());

#pragma omp critical
	{
		outfile.write(buffer, pos);
	}
}

std::string TrajectoryOutput::getDescription() const {
	return "Trajectory output";
}

ConditionalOutput::ConditionalOutput(std::string filename, std::string propName) :
		removeProperty(false) {
	propertyName = propName;
	outfile.open(filename.c_str());
	outfile
			<< "id, x, y, z, E, phi, theta, distance, i_id, i_x, i_y, i_z, i_E, i_phi, i_theta"
			<< std::endl;
}

ConditionalOutput::~ConditionalOutput() {
	outfile.close();
}

void ConditionalOutput::setRemoveProperty(bool removeProperty) {
	this->removeProperty = removeProperty;
}

void ConditionalOutput::process(Candidate *candidate) const {
	if (candidate->hasProperty(propertyName)) {
		char buffer[256];
		size_t pos = 0;

		pos += ::sprintf(buffer + pos, "%d", candidate->current.getId());
		Vector3 position = candidate->current.getPosition() / Mpc;
		pos += ::sprintf(buffer + pos, ", %f, %f, %f", position.x(),
				position.y(), position.z());
		const Vector3 &dir = candidate->current.getDirection();
		pos += ::sprintf(buffer + pos, ", %f, %f, %f",
				candidate->current.getEnergy(), dir.phi(), dir.theta());

		pos += ::sprintf(buffer + pos, ", %f",
				candidate->getTrajectoryLength());

		pos += ::sprintf(buffer + pos, ", %d", candidate->initial.getId());
		Vector3 ipos = candidate->initial.getPosition() / Mpc;
		pos += ::sprintf(buffer + pos, ", %f, %f, %f", ipos.x(), ipos.y(),
				ipos.z());
		const Vector3 &idir = candidate->initial.getDirection();
		pos += ::sprintf(buffer + pos, ", %f, %f, %f\n",
				candidate->initial.getEnergy(), idir.phi(), idir.theta());

#pragma omp critical
		{
			outfile.write(buffer, pos);
		}

		if (removeProperty) {
			candidate->removeProperty(propertyName);
		}
	}
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
		std::cout << candidate->current.getPosition() / Mpc << " Mpc, ";
		std::cout << candidate->current.getDirection().phi() << " / ";
		std::cout << candidate->current.getDirection().theta();
		std::cout << std::endl;
	}
}

std::string ShellOutput::getDescription() const {
	return "ShellOutput";
}

} // namespace mpc

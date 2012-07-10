#include "mpc/module/Output.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>

namespace mpc {

TrajectoryOutput::TrajectoryOutput(std::string name) {
	outfile.open(name.c_str());
	outfile << "# Age, HepId, E, posX, posY, posZ, dirX, dirY, dirZ, event\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	outfile.close();
}

void TrajectoryOutput::process(Candidate *candidate) const {
	char buffer[1024];
	size_t pos = 0;
	pos += ::sprintf(buffer + pos, "%.15g, %d, %.15g",
			candidate->getTrajectoryLength() / Mpc, candidate->current.getId(),
			candidate->current.getEnergy() / EeV);
	Vector3d position = candidate->current.getPosition() / Mpc;
	pos += ::sprintf(buffer + pos, ", %.15g, %.15g, %.15g", position.x, position.y,
			position.z);
	const Vector3d &dir = candidate->current.getDirection();
	pos += ::sprintf(buffer + pos, ", %.15g, %.15g, %.15g\n", dir.x, dir.y, dir.z);

#pragma omp critical
	{
		outfile.write(buffer, pos);
	}
}

std::string TrajectoryOutput::getDescription() const {
	return "Trajectory output";
}

ConditionalOutput::ConditionalOutput(std::string filename,
		std::string propName) {
	removeProperty = false;
	propertyName = propName;
	outfile.open(filename.c_str());
	outfile
			<< "# id, x, y, z, E, phi, theta, distance, i_id, i_x, i_y, i_z, i_E, i_phi, i_theta"
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
		size_t p = 0;

		p += ::sprintf(buffer + p, "%d", candidate->current.getId());

		const Vector3d &pos = candidate->current.getPosition() / Mpc;
		p += ::sprintf(buffer + p, ", %f, %f, %f", pos.x, pos.y, pos.z);

		const Vector3d &dir = candidate->current.getDirection();
		p += ::sprintf(buffer + p, ", %f, %f, %f",
				candidate->current.getEnergy() / EeV, dir.getPhi(),
				dir.getTheta());

		p += ::sprintf(buffer + p, ", %f",
				candidate->getTrajectoryLength() / Mpc);

		p += ::sprintf(buffer + p, ", %d", candidate->initial.getId());

		const Vector3d &ipos = candidate->initial.getPosition() / Mpc;
		p += ::sprintf(buffer + p, ", %f, %f, %f", ipos.x, ipos.y, ipos.z);

		const Vector3d &idir = candidate->initial.getDirection();
		p += ::sprintf(buffer + p, ", %f, %f, %f\n",
				candidate->initial.getEnergy() / EeV, idir.getPhi(),
				idir.getTheta());

#pragma omp critical
		{
			outfile.write(buffer, p);
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
		std::cout << candidate->current.getDirection().getPhi() << " / ";
		std::cout << candidate->current.getDirection().getTheta();
		std::cout << std::endl;
	}
}

std::string ShellOutput::getDescription() const {
	return "ShellOutput";
}

} // namespace mpc

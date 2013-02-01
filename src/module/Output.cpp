#include "mpc/module/Output.h"

#include <iomanip>
#include <stdio.h>

namespace mpc {

void ShellOutput::process(Candidate* c) const {
#pragma omp critical
	{
		std::cout << std::fixed << std::showpoint << std::setprecision(3)
				<< std::setw(6);
		std::cout << c->getTrajectoryLength() / Mpc << " Mpc,  ";
		std::cout << c->getRedshift() << ",  ";
		std::cout << c->current.getId() << ",  ";
		std::cout << c->current.getEnergy() / EeV << " EeV,  ";
		std::cout << c->current.getPosition() / Mpc << " Mpc,  ";
		std::cout << c->current.getDirection().getPhi() << " ";
		std::cout << c->current.getDirection().getTheta();
		std::cout << std::endl;
	}
}

std::string ShellOutput::getDescription() const {
	return "Shell output";
}

TrajectoryOutput::TrajectoryOutput(std::string name) {
	setDescription("Trajectory output");
	outfile.open(name.c_str());
	outfile << "# Age[Mpc]\t";
	outfile << "PDG_Code\t";
	outfile << "Energy[EeV]\t";
	outfile << "Position(X,Y,Z)[Mpc]\t";
	outfile << "Direction(X,Y,Z)\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	outfile.close();
}

void TrajectoryOutput::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%9.4f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	Vector3d pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%8.4f\t%8.4f\t%8.4f\t", pos.x, pos.y, pos.z);
	const Vector3d &dir = c->current.getDirection();
	p += sprintf(buffer + p, "%7.4f\t%7.4f\t%7.4f\n", dir.x, dir.y, dir.z);

#pragma omp critical
	{
		outfile.write(buffer, p);
		outfile.flush();
	}
}

ConditionalOutput::ConditionalOutput(std::string filename, std::string propName,
		bool removeProperty) :
		removeProperty(removeProperty), condition(propName) {

	setDescription(
			"ConditionalOutput, condition: " + propName + ", filename: "
					+ filename);

	outfile.open(filename.c_str());

	outfile << "# PDG_Code\t";
	outfile << "Energy[EeV]\t";
	outfile << "Position(X,Y,Z)[Mpc]\t";
	outfile << "Direction(Phi,Theta)\t";
	outfile << "Comoving distance [Mpc]\t";
	outfile << "Initial_PDG_Code\t";
	outfile << "Initial_Energy[EeV]\t";
	outfile << "Initial_Position(X,Y,Z)[Mpc]\t";
	outfile << "Initial_Direction(Phi,Theta)\n";
}

ConditionalOutput::~ConditionalOutput() {
	outfile.close();
}

void ConditionalOutput::setRemoveProperty(bool b) {
	removeProperty = b;
}

void ConditionalOutput::process(Candidate *c) const {
	if (not (c->hasProperty(condition)))
		return;

	if (removeProperty)
		c->removeProperty(condition);

	char buffer[256];
	size_t p = 0;

	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	const Vector3d &pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", pos.x, pos.y, pos.z);
	const Vector3d &dir = c->current.getDirection();
	p += sprintf(buffer + p, "%7.4f\t%7.4f\t", dir.getPhi(), dir.getTheta());
	p += sprintf(buffer + p, "%9.4f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->initial.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->initial.getEnergy() / EeV);
	const Vector3d &ipos = c->initial.getPosition() / Mpc;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", ipos.x, ipos.y, ipos.z);
	const Vector3d &idir = c->initial.getDirection();
	p += sprintf(buffer + p, "%7.4f\t%7.4f\n", idir.getPhi(), idir.getTheta());

#pragma omp critical
	{
		outfile.write(buffer, p);
		outfile.flush();
	}
}

TrajectoryOutput1D::TrajectoryOutput1D(std::string filename) {
	setDescription("Trajectory output");
	outfile.open(filename.c_str());
	outfile << "# Position(X)[Mpc]\t";
	outfile << "PDG_Code\t";
	outfile << "Energy[EeV]\n";
}

TrajectoryOutput1D::~TrajectoryOutput1D() {
	outfile.close();
}

void TrajectoryOutput1D::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;
	p += sprintf(buffer + p, "%8.4f\t", c->current.getPosition().x / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->current.getEnergy() / EeV);
#pragma omp critical
	{
		outfile.write(buffer, p);
		outfile.flush();
	}
}

EventOutput1D::EventOutput1D(std::string filename) {
	setDescription("Conditional output, Filename: " + filename);
	outfile.open(filename.c_str());
	outfile << "# PDG_Code\t";
	outfile << "Energy[EeV]\t";
	outfile << "Age[Mpc]";
	outfile << "Initial_PDG_Code\t";
	outfile << "Initial_Energy[EeV]\n";
}

EventOutput1D::~EventOutput1D() {
	outfile.close();
}

void EventOutput1D::process(Candidate *c) const {
	if (not (c->hasProperty("Detected")))
		return;

	c->removeProperty("Detected");

	char buffer[256];
	size_t p = 0;

	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%9.4f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->initial.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->initial.getEnergy() / EeV);

#pragma omp critical
	{
		outfile.write(buffer, p);
		outfile.flush();
	}
}

} // namespace mpc

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
	fout.open(name.c_str());
	fout << "# D\tID\tE\tX\tY\tZ\tPx\tPy\tPz\n";
	fout << "#\n";
	fout << "# D           Trajectory length\n";
	fout << "# ID          Particle type (PDG MC numbering scheme)\n";
	fout << "# E           Energy [EeV]\n";
	fout << "# X, Y, Z     Position [Mpc]\n";
	fout << "# Px, Py, Pz  Heading (unit vector of momentum)\n";
	fout << "#\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	fout.close();
}

void TrajectoryOutput::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.3f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	Vector3d pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%8.4f\t%8.4f\t%8.4f\t", pos.x, pos.y, pos.z);
	const Vector3d &dir = c->current.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\n", dir.x, dir.y, dir.z);

#pragma omp critical
	{
		fout.write(buffer, p);
		fout.flush();
	}
}

ConditionalOutput::ConditionalOutput(std::string fname, std::string cond) :
		condition(cond) {
	setDescription(
			"Conditional output, condition: " + cond + ", filename: " + fname);
	fout.open(fname.c_str());
	fout << "# D\tID\tID0\tE\tE0\tX\tY\tZ\tX0\tY0\tZ0\tPx\tPy\tPz\tP0x\tP0y\tP0z\n";
	fout << "#\n";
	fout << "# D           Trajectory length [Mpc]\n";
	fout << "# ID          Particle type (PDG MC numbering scheme)\n";
	fout << "# E           Energy [EeV]\n";
	fout << "# X, Y, Z     Position [Mpc]\n";
	fout << "# Px, Py, Pz  Heading (unit vector of momentum)\n";
	fout << "# Initial state: ID0, E0, ...\n";
	fout << "#\n";
}

ConditionalOutput::~ConditionalOutput() {
	fout.close();
}

void ConditionalOutput::process(Candidate *c) const {
	if (not (c->hasProperty(condition)))
		return;

	c->removeProperty(condition);

	char buffer[256];
	size_t p = 0;

	p += sprintf(buffer + p, "%8.3f\t", c->getTrajectoryLength() / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%10i\t", c->initial.getId());
	p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
	p += sprintf(buffer + p, "%8.4f\t", c->initial.getEnergy() / EeV);
	Vector3d pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", pos.x, pos.y, pos.z);
	Vector3d ipos = c->initial.getPosition() / Mpc;
	p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", ipos.x, ipos.y, ipos.z);
	Vector3d dir = c->current.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\t", dir.x, dir.y, dir.z);
	Vector3d idir = c->initial.getDirection();
	p += sprintf(buffer + p, "%8.5f\t%8.5f\t%8.5f\n", idir.x, idir.y, idir.z);

#pragma omp critical
	{
		fout.write(buffer, p);
		fout.flush();
	}
}

TrajectoryOutput1D::TrajectoryOutput1D(std::string filename) {
	setDescription("TrajectoryOutput, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#X\tID\tE\n";
	fout << "#\n";
	fout << "# X  Position [Mpc]\n";
	fout << "# ID Particle type\n";
	fout << "# E  Energy [EeV]\n";
}

TrajectoryOutput1D::~TrajectoryOutput1D() {
	fout.close();
}

void TrajectoryOutput1D::process(Candidate *c) const {
	char buffer[1024];
	size_t p = 0;
	p += sprintf(buffer + p, "%8.4f\t", c->current.getPosition().x / Mpc);
	p += sprintf(buffer + p, "%10i\t", c->current.getId());
	p += sprintf(buffer + p, "%8.4f\n", c->current.getEnergy() / EeV);
#pragma omp critical
	{
		fout.write(buffer, p);
		fout.flush();
	}
}

EventOutput1D::EventOutput1D(std::string filename) {
	setDescription("Conditional output, filename: " + filename);
	fout.open(filename.c_str());
	fout << "#ID\tE\tD\tID0\tE0\n";
	fout << "#\n";
	fout << "# ID  Particle type\n";
	fout << "# E   Energy [EeV]\n";
	fout << "# D   Comoving trajectory length [Mpc]\n";
	fout << "# ID0 Initial particle type\n";
	fout << "# E0  Initial energy [EeV]\n";
}

EventOutput1D::~EventOutput1D() {
	fout.close();
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
		fout.write(buffer, p);
		fout.flush();
	}
}

} // namespace mpc

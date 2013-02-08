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
	fout << "# D\tID\E\tX\tY\tZ\tPhi\tTheta\n";
	fout << "#\n";
	fout << "# D          Comoving trajectory length\n";
	fout << "# ID         Particle type\n";
	fout << "# E          Energy [EeV]\n";
	fout << "# X, Y, Z    Position [Mpc]\n";
	fout << "# Phi, Theta Direction\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	fout.close();
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
	p += sprintf(buffer + p, "%7.4f\t%7.4f\n", dir.getPhi(), dir.getTheta());

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
	fout << "# ID\tE\tY\tY\tZ\tPhi\tTheta\tD\tID0\tE0\tX0\tY0\tZ0\tPhi0\tTheta0\n";
	fout << "#\n";
	fout << "# Final state:\n";
	fout << "#    ID           Particle type\n";
	fout << "#    E            Energy [EeV]\n";
	fout << "#    X, Y, Z      Position [Mpc]\n";
	fout << "#    Phi, Theta   Direction\n";
	fout << "#    D            Comoving trajectory length [Mpc]\n";
	fout << "#\n";
	fout << "# Initial state:\n";
	fout << "#    I0           Particle type\n";
	fout << "#    E0           Energy [EeV]\n";
	fout << "#    X0, Y0, Z0   Position [Mpc]\n";
	fout << "#    Phi0, Theta0 Direction\n";
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

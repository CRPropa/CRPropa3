#include "mpc/module/Output.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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
		removeProperty(removeProperty) {
	setDescription(
			"Conditional output, condition: " + propName + ", filename: "
					+ filename);
	condition = propName;
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
	if (c->hasProperty(condition)) {
		if (removeProperty)
			c->removeProperty(condition);

		char buffer[256];
		size_t p = 0;

		p += sprintf(buffer + p, "%10i\t", c->current.getId());
		p += sprintf(buffer + p, "%8.4f\t", c->current.getEnergy() / EeV);
		const Vector3d &pos = c->current.getPosition() / Mpc;
		p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", pos.x, pos.y, pos.z);
		const Vector3d &dir = c->current.getDirection();
		p += sprintf(buffer + p, "%7.4f\t%7.4f\t", dir.getPhi(),
				dir.getTheta());
		p += sprintf(buffer + p, "%9.4f\t", c->getTrajectoryLength() / Mpc);
		p += sprintf(buffer + p, "%10i\t", c->initial.getId());
		p += sprintf(buffer + p, "%8.4f\t", c->initial.getEnergy() / EeV);
		const Vector3d &ipos = c->initial.getPosition() / Mpc;
		p += sprintf(buffer + p, "%9.4f\t%9.4f\t%9.4f\t", ipos.x, ipos.y,
				ipos.z);
		const Vector3d &idir = c->initial.getDirection();
		p += sprintf(buffer + p, "%7.4f\t%7.4f\n", idir.getPhi(),
				idir.getTheta());

#pragma omp critical
		{
			outfile.write(buffer, p);
			outfile.flush();
		}
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
	if (c->isActive())
		return;
	if (c->hasProperty("Detected")) {
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
}

CRPropa2EventOutput::CRPropa2EventOutput(std::string filename) {
	setDescription("Event output in CRPropa2 format");
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file\n";
	outfile
			<< "#Format - Particle_Type Initial_Particle_Type Initial_Position[X,Y,Z](Mpc) Initial_Momentum[E,theta,phi](EeV) Time(Mpc) Position[X,Y,Z](Mpc) Momentum[E,theta,phi](EeV)\n";
}

CRPropa2EventOutput::~CRPropa2EventOutput() {
	outfile.close();
}

void CRPropa2EventOutput::process(Candidate *candidate) const {
	if (candidate->hasProperty("Detected")) {
		char buffer[256]; // max. 256 characters per line
		size_t p = 0; // length of line

		int cid = convertToCRPropaId(candidate->current.getId());
		p += sprintf(buffer + p, "%i ", cid);

		int iid = convertToCRPropaId(candidate->initial.getId());
		p += sprintf(buffer + p, "%i ", iid);

		const Vector3d &ipos = candidate->initial.getPosition() / Mpc;
		p += sprintf(buffer + p, "%.4f %.4f %.4f ", ipos.x, ipos.y, ipos.z);

		double iPhi = candidate->initial.getDirection().getPhi();
		double iTheta = candidate->initial.getDirection().getTheta();
		double iE = candidate->initial.getEnergy() / EeV;
		p += sprintf(buffer + p, "%.4f %.4f %.4f ", iE, iPhi, iTheta);

		double t = candidate->getTrajectoryLength() / Mpc;
		p += sprintf(buffer + p, "%.4f ", t);

		const Vector3d &pos = candidate->current.getPosition() / Mpc;
		p += sprintf(buffer + p, "%.4f %.4f %.4f ", pos.x, pos.y, pos.z);

		double phi = candidate->current.getDirection().getPhi();
		double theta = candidate->current.getDirection().getTheta();
		double E = candidate->current.getEnergy() / EeV;
		p += sprintf(buffer + p, "%.4f %.4f %.4f\n", E, phi, theta);

#pragma omp critical
		{
			outfile.write(buffer, p);
			outfile.flush();
		}
	}
}

CRPropa2TrajectoryOutput::CRPropa2TrajectoryOutput(std::string filename) {
	setDescription("Trajectory output in CRPropa2 format");
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file" << std::endl
			<< "#Format - Particle_Type Initial_Particle_Type Time(Mpc) Position[X,Y,Z](Mpc) Momentum[X,Y,Z](EeV) Energy(EeV)"
			<< std::endl;
}

CRPropa2TrajectoryOutput::~CRPropa2TrajectoryOutput() {
	outfile.close();
}

void CRPropa2TrajectoryOutput::process(Candidate *candidate) const {
	char buffer[256];
	size_t p = 0;

	int cid = convertToCRPropaId(candidate->current.getId());
	p += sprintf(buffer + p, "%i ", cid);

	int iid = convertToCRPropaId(candidate->initial.getId());
	p += sprintf(buffer + p, "%i ", iid);

	double t = candidate->getTrajectoryLength() / Mpc;
	p += sprintf(buffer + p, "%.4f ", t);

	const Vector3d &pos = candidate->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%.4f %.4f %.4f ", pos.x, pos.y, pos.z);

	const Vector3d &mom = candidate->current.getMomentum() / EeV;
	p += sprintf(buffer + p, "%.4g %.4g %.4g ", mom.x, mom.y, mom.z);

	double E = candidate->current.getEnergy() / EeV;
	p += sprintf(buffer + p, "%.4f\n", E);

#pragma omp critical
	{
		outfile.write(buffer, p);
		outfile.flush();
	}
}

CRPropa2TrajectoryOutput1D::CRPropa2TrajectoryOutput1D(std::string filename) {
	setDescription("Trajectory output");
	outfile.open(filename.c_str());
	outfile << "# Position(X)[Mpc]\t";
	outfile << "PDG_Code\t";
	outfile << "Energy[EeV]\n";
}

CRPropa2TrajectoryOutput1D::~CRPropa2TrajectoryOutput1D() {
	outfile.close();
}

void CRPropa2TrajectoryOutput1D::process(Candidate *c) const {
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

CRPropa2EventOutput1D::CRPropa2EventOutput1D(std::string filename) {
	setDescription("Conditional output, Filename: " + filename);
	outfile.open(filename.c_str());
	outfile << "# PDG_Code\t";
	outfile << "Energy[EeV]\t";
	outfile << "Age[Mpc]";
	outfile << "Initial_PDG_Code\t";
	outfile << "Initial_Energy[EeV]\n";
}

CRPropa2EventOutput1D::~CRPropa2EventOutput1D() {
	outfile.close();
}

void CRPropa2EventOutput1D::process(Candidate *c) const {
	if (c->isActive())
		return;
	if (c->hasProperty("Detected")) {
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
}

} // namespace mpc

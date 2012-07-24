#include "mpc/module/Output.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>

namespace mpc {

TrajectoryOutput::TrajectoryOutput(std::string name) {
	setDescription("Trajectory output");
	outfile.open(name.c_str());
	outfile << "# Age, HepId, E, posX, posY, posZ, dirX, dirY, dirZ, event\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	outfile.close();
}

void TrajectoryOutput::process(Candidate *candidate) const {
	char buffer[1024];
	size_t pos = 0;
	pos += ::sprintf(buffer + pos, "%.4f, %d, %.4f",
			candidate->getTrajectoryLength() / Mpc, candidate->current.getId(),
			candidate->current.getEnergy() / EeV);
	Vector3d position = candidate->current.getPosition() / Mpc;
	pos += ::sprintf(buffer + pos, ", %.4f, %.4f, %.4f", position.x, position.y,
			position.z);
	const Vector3d &dir = candidate->current.getDirection();
	pos += ::sprintf(buffer + pos, ", %.4f, %.4f, %.4f\n", dir.x, dir.y, dir.z);

#pragma omp critical
	{
		outfile.write(buffer, pos);
	}
}

ConditionalOutput::ConditionalOutput(std::string filename,
		std::string propName) {
	setDescription("ConditionalOutput, condition: " + propName);
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

		p += sprintf(buffer + p, "%i", candidate->current.getId());

		const Vector3d &pos = candidate->current.getPosition() / Mpc;
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f", pos.x, pos.y, pos.z);

		const Vector3d &dir = candidate->current.getDirection();
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f",
				candidate->current.getEnergy() / EeV, dir.getPhi(),
				dir.getTheta());

		p += sprintf(buffer + p, ", %.4f",
				candidate->getTrajectoryLength() / Mpc);

		p += sprintf(buffer + p, ", %i", candidate->initial.getId());

		const Vector3d &ipos = candidate->initial.getPosition() / Mpc;
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f", ipos.x, ipos.y, ipos.z);

		const Vector3d &idir = candidate->initial.getDirection();
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f\n",
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

//  _fDataStream << "#CRPropa - Output data file" << endl ;
//  if (lUnivType == "One Dimension" && lRecordType == "Events") {
//    _fDataStream << "#Format - Particle_Type Initial_Particle_Type Initial_Position(Mpc) Initial_Energy(EeV) Time(Mpc) Energy(EeV)" << endl ;
//  } else if (lUnivType == "One Dimension" && lRecordType == "Full Trajectories") {
//    _fDataStream << "#Format - Particle_Type Initial_Particle_Type Time(Mpc) Position(Mpc) Energy(EeV)" << endl ;
//  } else throw TCrpErr("Output format determination failed") ;

CRPropa2EventOutput::CRPropa2EventOutput(std::string filename) {
	setDescription("Event output in CRPropa2 format");
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file" << std::endl
			<< "#Format - Particle_Type Initial_Particle_Type Initial_Position[X,Y,Z](Mpc) Initial_Momentum[E,theta,phi](EeV) Time(Mpc) Position[X,Y,Z](Mpc) Momentum[E,theta,phi](EeV)"
			<< std::endl;
}

CRPropa2EventOutput::~CRPropa2EventOutput() {
	outfile.close();
}

void CRPropa2EventOutput::process(Candidate *candidate) const {
	if (candidate->isActive())
		return;
	if (candidate->hasProperty("Detected")) {
		char buffer[256]; // max. 256 characters per line
		size_t p = 0; // length of line

		int cid = convertToCRPropaId(candidate->current.getId());
		p += sprintf(buffer + p, "%i, ", cid);

		int iid = convertToCRPropaId(candidate->initial.getId());
		p += sprintf(buffer + p, "%i, ", iid);

		const Vector3d &ipos = candidate->initial.getPosition() / Mpc;
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f", ipos.x, ipos.y, ipos.z);

		double iPhi = candidate->initial.getDirection().getPhi();
		double iTheta = candidate->initial.getDirection().getTheta();
		double iE = candidate->initial.getEnergy() / EeV;
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f", iE, iPhi, iTheta);

		double t = candidate->getTrajectoryLength() / Mpc;
		p += sprintf(buffer + p, ", %.4f", t);

		const Vector3d &pos = candidate->current.getPosition() / Mpc;
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f", ipos.x, ipos.y, ipos.z);

		double phi = candidate->current.getDirection().getPhi();
		double theta = candidate->current.getDirection().getTheta();
		double E = candidate->current.getEnergy() / EeV;
		p += sprintf(buffer + p, ", %.4f, %.4f, %.4f\n", E, phi, theta);

#pragma omp critical
		{
			outfile.write(buffer, p);
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
	p += sprintf(buffer + p, "%i", cid);

	int iid = convertToCRPropaId(candidate->initial.getId());
	p += sprintf(buffer + p, ", %i", iid);

	double t = candidate->getTrajectoryLength() / Mpc;
	p += sprintf(buffer + p, ", %.4f", t);

	const Vector3d &pos = candidate->current.getPosition() / Mpc;
	p += sprintf(buffer + p, ", %.4f, %.4f, %.4f", pos.x, pos.y, pos.z);

	const Vector3d &mom = candidate->current.getMomentum() / EeV;
	p += sprintf(buffer + p, ", %.4g, %.4g, %.4g", mom.x, mom.y, mom.z);

	double E = candidate->current.getEnergy() / EeV;
	p += sprintf(buffer + p, ", %.4f\n", E);

#pragma omp critical
	{
		outfile.write(buffer, p);
	}
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
	return "Shell output";
}

} // namespace mpc

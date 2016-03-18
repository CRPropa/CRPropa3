#include "crpropa/module/OutputCRPropa2.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Cosmology.h"

namespace crpropa {

CRPropa2EventOutput3D::CRPropa2EventOutput3D(std::string filename) {
	setDescription("Event output in CRPropa2 format");
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file\n"
			<< "#Format - Particle_Type "
			<< "Initial_Particle_Type "
			<< "Initial_Position[X,Y,Z](Mpc) "
			<< "Initial_Momentum[E(EeV),theta,phi] "
			<< "Time(Mpc, light travel distance) "
			<< "Position[X,Y,Z](Mpc) "
			<< "Momentum[E(EeV),theta,phi]\n";
}

CRPropa2EventOutput3D::~CRPropa2EventOutput3D() {
	outfile.close();
}

void CRPropa2EventOutput3D::process(Candidate *c) const {
	char buffer[256]; // max. 256 characters per line
	size_t p = 0; // length of line

	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->current.getId()));
	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->source.getId()));

	const Vector3d &ipos = c->source.getPosition() / Mpc;
	p += sprintf(buffer + p, "%.4f %.4f %.4f ", ipos.x, ipos.y, ipos.z);

	double iPhi = c->source.getDirection().getPhi();
	double iTheta = c->source.getDirection().getTheta();
	double iE = c->source.getEnergy() / EeV;
	p += sprintf(buffer + p, "%.4f %.4f %.4f ", iE, iPhi, iTheta);

	double t = comoving2LightTravelDistance(c->getTrajectoryLength()) / Mpc;
	p += sprintf(buffer + p, "%.4f ", t);

	const Vector3d &pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%.4f %.4f %.4f ", pos.x, pos.y, pos.z);

	double phi = c->current.getDirection().getPhi();
	double theta = c->current.getDirection().getTheta();
	double E = c->current.getEnergy() / EeV;
	p += sprintf(buffer + p, "%.4f %.4f %.4f\n", E, phi, theta);

#pragma omp critical
	outfile.write(buffer, p);
}

void CRPropa2EventOutput3D::close() {
	outfile.flush();
}

CRPropa2TrajectoryOutput3D::CRPropa2TrajectoryOutput3D(std::string filename) {
	setDescription("Trajectory output in CRPropa2 format");
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file\n"
			<< "#Format - Particle_Type "
			<< "Initial_Particle_Type "
			<< "Time(Mpc, light travel distance) "
			<< "Position[X,Y,Z](Mpc) "
			<< "Momentum[X(EeV),Y,Z] "
			<< "Energy(EeV)\n";
}

CRPropa2TrajectoryOutput3D::~CRPropa2TrajectoryOutput3D() {
	outfile.close();
}

void CRPropa2TrajectoryOutput3D::process(Candidate *c) const {
	char buffer[256];
	size_t p = 0;

	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->current.getId()));
	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->source.getId()));

	double t = comoving2LightTravelDistance(c->getTrajectoryLength()) / Mpc;
	p += sprintf(buffer + p, "%.4f ", t);

	const Vector3d &pos = c->current.getPosition() / Mpc;
	p += sprintf(buffer + p, "%.4f %.4f %.4f ", pos.x, pos.y, pos.z);

	const Vector3d &mom = c->current.getMomentum() / EeV;
	p += sprintf(buffer + p, "%.4g %.4g %.4g ", mom.x, mom.y, mom.z);

	p += sprintf(buffer + p, "%.4f\n", c->current.getEnergy() / EeV);

#pragma omp critical
	outfile.write(buffer, p);
}

void CRPropa2TrajectoryOutput3D::close() {
	outfile.flush();
}

CRPropa2TrajectoryOutput1D::CRPropa2TrajectoryOutput1D(std::string filename) {
	setDescription("Trajectory output");
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file\n"
			<< "#Format - Position(Mpc) "
			<< "Particle_Type "
			<< "Energy(EeV)\n";
}

CRPropa2TrajectoryOutput1D::~CRPropa2TrajectoryOutput1D() {
	outfile.close();
}

void CRPropa2TrajectoryOutput1D::process(Candidate *c) const {
	char buffer[256];
	size_t p = 0;

	p += sprintf(buffer + p, "%.4f ", comoving2LightTravelDistance(c->current.getPosition().x) / Mpc);
	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->current.getId()));
	p += sprintf(buffer + p, "%.4f\n", c->current.getEnergy() / EeV);

#pragma omp critical
	outfile.write(buffer, p);
}

void CRPropa2TrajectoryOutput1D::close() {
	outfile.flush();
}

CRPropa2EventOutput1D::CRPropa2EventOutput1D(std::string filename) {
	setDescription("Conditional output, Filename: " + filename);
	outfile.open(filename.c_str());
	outfile << "#CRPropa - Output data file\n"
			<< "#Format - Energy(EeV) "
			<< "Time(Mpc, light travel distance) "
			<< "Initial_Particle_Type "
			<< "Initial_Energy(EeV)\n";
}

CRPropa2EventOutput1D::~CRPropa2EventOutput1D() {
	outfile.close();
}

void CRPropa2EventOutput1D::process(Candidate *c) const {
	char buffer[256];
	size_t p = 0;

	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->current.getId()));
	p += sprintf(buffer + p, "%.4f ", c->current.getEnergy() / EeV);
	double t = comoving2LightTravelDistance(c->getTrajectoryLength()) / Mpc;
	p += sprintf(buffer + p, "%.4f ", t);
	p += sprintf(buffer + p, "%i ", convertToCRPropa2NucleusId(c->source.getId()));
	p += sprintf(buffer + p, "%.4f\n", c->source.getEnergy() / EeV);

#pragma omp critical
	outfile.write(buffer, p);
}

void CRPropa2EventOutput1D::close() {
	outfile.flush();
}

} // namespace crpropa
